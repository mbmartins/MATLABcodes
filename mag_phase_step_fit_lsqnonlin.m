% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

%signal generation
F0 = 60; %Hz
F1 = 60; %Hz
SampleRate = 4800; %Hz
dt = 1/SampleRate;
AnalysisCycles = 6;
NSamples = floor(AnalysisCycles*SampleRate/F0);
n = -NSamples/2:(NSamples/2-1); %discrete time vector

Vm = 100; %70*sqrt(2);
Xm = Vm;
Ps = 120; %phase in degrees
% Phase in radians
Ph = Ps*pi/180;
KaS = 0;   % IEEE Std phase (angle) step index: 10 degrees
KxS = 0.1;   % magnitude step index: 0.1 
Wf = 2*pi*F1;  % fundamental frequency
SNR = 90.5; %dB SNR = 20 log_10 Asinal/Aruido => Aruido = Asinal/10^(SNR/20)
Aruido = Vm/10^(SNR/20);
    
for ti = 1:9
    
    tau_pp = 0.1*ti; % relative time of step in percent of total time 
    tau_0 = (tau_pp - 0.5)*NSamples; %discrete time displacement
    n = -NSamples/2:(NSamples/2-1); %discrete time vector
    n = n - tau_0;
    t = n*dt; %time vector

    %%%% Signal parameters estimation
    
    u = zeros(length(Xm),length(t));
    tau = 0;
    u(length(Xm),t >= tau) = u(length(Xm),t >= tau) + 1;
    
    %modelo 1
    % Phase step
    % Model: f(x) = x1*cos(w*t + phi + x2*u(t - tau)) 
    if KaS ~= 0
        f = @(x) x(1)*cos(x(2)*t + x(3) + x(4)*(pi/180)*u);
        xnom1 = [Vm 2*pi*F1 Ph KaS];
    else
        %mag
        %f = @(x) x(1)*(1+x(4)*u).*cos(x(2)*t + x(3));
        %xnom1 = [Vm 2*pi*F1 Ph KxS];
        f = @(x) x(1)*(1+x(2)*u).*cos(x(3)*t + x(4));
        xnom1 = [Vm KxS 2*pi*F1 Ph];
    end
    
    err = @(x) (Signal - f(x)).^2;

    % Nonlinear fit
    % Monte Carlo analysis
    xnom = xnom1;
    Niter = 1000;
    x0 = xnom;  %x0 is fixed - first guess are the nominal values
    
    for k = 1:Niter
        %first guess
        if KaS ~= 0
        % phase         X1  w    ph   x3 (KaS)
            par_var = [0.02 0.01 0.01 0.02]; % parameter variation in percent related to nominal
        else
        % mag          x1  x2(KxS)  wf    ph  
            par_var = [0.02 0.02     0.01  0.01]; % parameter variation in percent related to nominal            
        end

        rng('shuffle');
        rn = (rand(1,length(xnom1))-0.5);
        xr = xnom1.*(1+(par_var/100).*rn);
        if KaS ~= 0
            if xr(3) == 0   %correção para fase zero
                xr(3) = xr(3) + (par_var(3)/100)*rn(3)*120; 
            end
        end
        utau = 2;  %number of dts - incerteza da estimação de tau 
        u = zeros(length(Xm),length(t));
        tau = dt*randi([-utau utau],1);
        u(length(Xm),t >= tau) = u(length(Xm),t >= tau) + 1;
        
        %regenerates the signal with the uncertainties
            i = 1;
            Xm = xr(1); 
            if KaS ~= 0
                Wf(i) = xr(2); 
                Ph(i) = xr(3);
                KxS_(i) = 0;
                KaS_(i) = xr(4);
            else
                Ph(i) = xr(4); 
                Wf(i) = xr(3);
                KxS_(i) = xr(2);
                KaS_(i) = 0;
            end
            
            Ain = zeros(length(Xm),length(t));
            % Amplitude Step: applied after time passes 0
            Ain(i,:) = Xm(i);
            Ain(i,t >= 0) = Ain(i,t >= 0) * (1 + KxS_(i));
            %Phase step
            Theta(i,:) = (Wf(i)*t) ...                         % Fundamental
                             + Ph(i);               % phase shift
            Theta(i,t >= 0) = Theta(i,t >= 0) + (KaS_(i) * pi/180);
            cSignal = (Ain.*exp(-1i.*Theta));
            Signal = real(cSignal) + Aruido*(rand(1,length(t))-0.5);
            err = @(x) (Signal - f(x)).^2;
        y0 = f(x0);
        x = xnom1;
        %estimated signal
%         y = f(xnom1);
%         plot(t,Signal,'.b',t,y,'r')
     
        % figure
        % plot(t,Signal - y0,'r')
        % title('Point by point estimation error: Signal - y')

        % Levenberg-marquardt from optimization toolbox
        tol = 1e-7;
        OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','OptimalityTolerance',tol);
        OPTIONS.StepTolerance = 1e-12;
        [X,RESNORM,RESIDUAL,exitflag,output] = lsqnonlin(err,x0,[],[],OPTIONS);
        Y = f(X);
        errors(k,:) = (xr - X)./xr;
        if xr(3) == 0
           errors(k,:) = (xr - X)./120; 
        end
    end

    amp_rand(ti,:) = max(abs(rn));
    ERR_MAX = max(errors)*100';   %erros maximos em %
    ERR_MIN = min(errors)*100';   %erros minimos em %

    % figure
    % for k = 1:2
    %     subplot(2,2,k); hold on;
    %     histogram(errors(:,k));
    % end

    MEAN_ERR = mean(errors)*100;  %erros medios em %
    STDEV_ERR = std(errors)*100;  % desvio padrao em %
    %hold on; plot(X(2),X(1),'rx', x0(2), x0(1), 'ro')

    %Valores de magnitude X1
    Mag_errmax(ti,1) = ERR_MAX(1);
    Mag_errmin(ti,1) = ERR_MIN(1);
    Mag_errmed(ti,1) = MEAN_ERR(1);
    Mag_stddev(ti,1) = STDEV_ERR(1);
    
    %Valores de frequencia
    if KaS ~= 0
        ifreq = 2;
    else
        ifreq = 3;
    end
    Freq_errmax(ti,1) = ERR_MAX(ifreq);
    Freq_errmin(ti,1) = ERR_MIN(ifreq);
    Freq_errmed(ti,1) = MEAN_ERR(ifreq);
    Freq_stddev(ti,1) = STDEV_ERR(ifreq);
    
    %Valores de fase
    if KaS ~= 0
        iph = 3;
    else
        iph = 4;
    end
    Ph_errmax(ti,1) = ERR_MAX(iph);
    Ph_errmin(ti,1) = ERR_MIN(iph);
    Ph_errmed(ti,1) = MEAN_ERR(iph);
    Ph_stddev(ti,1) = STDEV_ERR(iph);    
    
    %Valores de degrau (KxS ou KaS)
    if KaS ~= 0
        ik = 4;
    else
        ik = 2;
    end
    K_errmax(ti,1) = ERR_MAX(ik);
    K_errmin(ti,1) = ERR_MIN(ik);
    K_errmed(ti,1) = MEAN_ERR(ik);
    K_stddev(ti,1) = STDEV_ERR(ik);    
    
%    Freq_err_per_std(ti,1) = ERR_MAX(2)/STDEV_ERR(2);
    
end    