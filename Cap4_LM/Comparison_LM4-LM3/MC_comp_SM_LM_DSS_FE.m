% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

TVE = @(Ve_r,Ve) 100*sqrt(((real(Ve_r) - real(Ve))^2 + (imag(Ve_r) - imag(Ve))^2 )/(real(Ve_r)^2 + imag(Ve_r)^2)) ;

%signal generation
F0 = 50; %Hz
F1 = 50; %Hz
SNR = 93.5; %dB
fs = 5000; %Hz
dt = 1/fs;
AnalysisCycles = 3; 
NSamples = floor(AnalysisCycles*fs/F0);
n = -NSamples/2:(NSamples/2-1); %discrete time vector
tau_0 = 0; %discrete time displacement
n = n + tau_0;
t = n*dt; %time vector
Vm = 9; %70*sqrt(2);
Ps = 360; %phase in degrees
% Phase in radians
Ph = Ps*pi/180;

Xref = (Vm/sqrt(2))*exp(-1i*Ph);
    

ha = 0;   % IEEE Std phase (angle) step index: 10 degrees
hm = 0.1;   % magnitude step index: 0.1 
Wf = 2*pi*F1;  % fundamental frequency

Xm = Vm; %for now, single phase; TODO: 6-channels
Ain = zeros(length(Xm),length(t));
% Amplitude Step: applied after time passes 0
i = 1;
Ain(i,:) = Xm(i);
Ain(i,t >= 0) = Ain(i,t >= 0) * (1 + hm(i));
%Phase step
Theta(i,:) = (Wf(i)*t) ...                         % Fundamental
                 + Ph(i);               % phase shift
Theta(i,t >= 0) = Theta(i,t >= 0) + (ha(i) * pi/180);
cSignal = (Ain.*exp(-1i.*Theta));
Signal = real(cSignal);

%%%% Signal parameters estimation
% Phase step
% Model: f(x) = x1*cos(w*t + phi + x2*u(t - tau))
u = zeros(length(Xm),length(t));
tau = 0;
u(length(Xm),t >= tau) = u(length(Xm),t >= tau) + 1;

        %modelo 1
        % Phase step
        % Model: f(x) = x1*cos(w*t + phi + x2*u(t - tau)) 
        if ha ~= 0
            %f = @(x) x(1)*cos(x(2)*t + x(3) + x(4)*(pi/180)*u);
            f = @(x) x(1)*cos(x(2)*t + x(3) + 2*pi + x(4)*(pi/180)*u);
            xnom = [Xm 2*pi*F1 Ph ha];
            par_var = [1 0.05 1 1]; % parameter variation in percent related to nominal
        else
            %mag
            %f = @(x) x(1)*(1+x(4)*u).*cos(x(2)*t + x(3));
            %xnom1 = [Vm 2*pi*F1 Ph KxS];
            f = @(x) x(1)*(1+x(2)*u).*cos(x(3)*t + x(4));
            xnom = [Xm hm 2*pi*F1 Ph];
            %par_var = [1 1 0.05 1]; % parameter variation in percent related to nominal
            par_var = [0 0.05 0 0];
        end

err = @(x) (Signal - f(x));
%first estimates for x
% SNR = 20 log_10 Asinal/Aruido => Aruido = Asinal/10^(SNR/20)
Aruido = Vm/10^(SNR/20);

%acertar x com o modelo escolhido
x = xnom;

% Nonlinear fit
% Monte Carlo analysis
Niter = 1000
x0 = xnom;  %x0 is fixed - first guess are the nominal values
for n = 1:Niter;
    
    xr = xnom.*(1+(par_var/100).*(rand(1,length(xnom))-0.5));
    
    %regenerates the signal with the uncertainties
        i = 1;
        if ha ~= 0 %modelo salto fase
            Xm = xr(1); Wf(i) = Wf; Ph(i) = xr(3); ha(i) = xr(4);   
        else
            Xm = xr(1); Wf(i) = xr(3); Ph(i) = xr(4); hm(i) = xr(2);
        end
        Xref = (Xm/sqrt(2))*exp(-1i*Ph);
            
        Ain = zeros(length(Xm),length(t));
        % Amplitude Step: applied after time passes 0
        Ain(i,:) = Xm(i);
        Ain(i,t >= 0) = Ain(i,t >= 0) * (1 + hm(i));
        %Phase step
        Theta(i,:) = (Wf(i)*t) ...                         % Fundamental
                         + Ph(i);               % phase shift
        Theta(i,t >= 0) = Theta(i,t >= 0) + (ha(i) * pi/180);
        cSignal = (Ain.*exp(-1i.*Theta));
        Signal = real(cSignal) + Aruido*(rand(1,length(t))-0.5);
        
        err = @(x) (Signal - f(x));
    y0 = f(x0);
    
    % Levenberg-marquardt from optimization toolbox
    tol = 1e-7;
    OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','OptimalityTolerance',tol);
    OPTIONS.StepTolerance = 1e-12;
    OPTIONS.Display = 'off';
    [X,RESNORM,RESIDUAL,exitflag,output,lambda,jacobian] = lsqnonlin(err,x0,[],[],OPTIONS);
    Y = f(X);
    
    % estimação por DSS
    [DSS_F_Hz,DSS_Mag(n), DSS_Fase_deg(n)] = DSS(F1, F0,AnalysisCycles, fs, Signal);
    DSS_FE(n) = DSS_F_Hz - F1;
    DSS_FE_max(n) = max(abs(DSS_FE));
    
    if ha ~= 0
        ferrHz(n) = X(2)/(2*pi) - F1;
    else
        ferrHz(n) = X(3)/(2*pi) - F1;
    end
    errors(n,:) = (xr - X)./xr;
    HLM4_FEmax(n) = max(abs(ferrHz));
    
    %TVE
    if ha ~= 0 % salto fase
        X_HLM4 = (X(1)/sqrt(2))*exp(-1i*X(3));
    else
        X_HLM4 = (X(1)/sqrt(2))*exp(-1i*X(4));
    end
    X_DSS = DSS_Mag(n)*exp(-1i*DSS_Fase_deg(n)*pi/180);
    
    HLM4_TVE(n) = TVE(Xref,X_HLM4);
    DSS_TVE(n) = TVE(Xref,X_DSS);
    
    HLM4_TVEmax(n) = max(HLM4_TVE);
    DSS_TVEmax(n) = max(DSS_TVE);
    
    fprintf("MC= "+ n + "; FEmax (HLM4= " + HLM4_FEmax(n) + "; DSS= " + DSS_FE_max(n) + ") Hz; \n");
end

ERR_MAX = max(abs(errors))*100';   %erros maximos em %


% 
% figure
% for k = 1:2
%     subplot(2,2,k); hold on;
%     histogram(errors(:,k));
% end
% 
% MEAN_ERR = mean(errors)
% STDEV_ERR = std(errors)
% %hold on; plot(X(2),X(1),'rx', x0(2), x0(1), 'ro')
% 
% %freq_err = ERR_MAX(2)*F1/(2*pi)
% ferrHz_max = max(abs(ferrHz))

% save('C:\Users\mbrit\Documents\MATLAB\MATLABcodes\Cap4_LM\Comparison_LM4-LM3\comp_LM_DSS')
    