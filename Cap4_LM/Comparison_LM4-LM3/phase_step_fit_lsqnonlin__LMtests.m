% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

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
Vm = 1; %70*sqrt(2);
Ps = 360; %phase in degrees

% Phase in radians
Ph = Ps*pi/180;

KaS = 10;   % IEEE Std phase (angle) step index: 10 degrees
KxS = 0.0;   % magnitude step index: 0.1 
Wf = 2*pi*F1;  % fundamental frequency

Xm = Vm; %for now, single phase; TODO: 6-channels
Ain = zeros(length(Xm),length(t));
% Amplitude Step: applied after time passes 0
i = 1;
Ain(i,:) = Xm(i);
Ain(i,t >= 0) = Ain(i,t >= 0) * (1 + KxS(i));
%Phase step
Theta(i,:) = (Wf(i)*t) ...                         % Fundamental
                 + Ph(i);               % phase shift
Theta(i,t >= 0) = Theta(i,t >= 0) + (KaS(i) * pi/180);
cSignal = (Ain.*exp(-1i.*Theta));
Signal = real(cSignal);

%%%% Signal parameters estimation
% Phase step
% Model: f(x) = x1*cos(w*t + phi + x2*u(t - tau))
u = zeros(length(Xm),length(t));
tau = 0;
u(length(Xm),t >= tau) = u(length(Xm),t >= tau) + 1;

%modelo 1 - LM4
f = @(x) x(1)*cos(x(2)*t + x(3) + x(4)*(pi/180)*u);
%modelo 2 - LM3 - frequencia fixa
% f = @(x) x(1)*cos(Wf*t + x(2) + x(3)*(pi/180)*u);
%modelo 3 - LM2 - frequencia e fase fixas
%f = @(x) x(1)*cos(Wf*t + Ph + x(2)*(pi/180)*u);


err = @(x) (Signal - f(x)).^2;
%first estimates for x
% SNR = 20 log_10 Asinal/Aruido => Aruido = Asinal/10^(SNR/20)
Aruido = Vm/10^(SNR/20);

xnom1 = [Vm 2*pi*F1 Ph KaS];
xnom2 = [Vm Ph KaS];
xnom3 = [Vm KaS];

%acertar x com o modelo escolhido
x = xnom1;

%estimated signal
y = f(x);
plot(t,Signal,'.b',t,y,'r')
% 
% %sum of squared errors with variation of x
% dx1 = 1e-2;
% x1 = (x(1)-50*dx1:dx1:x(1)+49*dx1);
% x2 = (x(2)-50*dx1:dx1:x(2)+49*dx1);
% 
% for i = 1:length(x1)   
% %    for j = 1:length(tau)
% for j = 1:length(x2)
%         %u = zeros(length(Xm),length(t));
%         %u(length(Xm),t >= tau(j)) = u(length(Xm),t >= tau(j)) + 1;
%         y = f([x1(i) x2(j)]) + Aruido*(rand(1,length(t))-0.5);
%         errq(i,j) = sum((Signal - y).^2);
% end
% %    end
% end
% 
% figure
% %plot(A,err)
% ES = surf(x2,x1,errq);
% title('Error surface: x2 x x1 x error')

% Nonlinear fit
% Monte Carlo analysis
xnom = xnom1;
Niter = 1000
x0 = xnom;  %x0 is fixed - first guess are the nominal values
for n = 1:Niter;
    %first guess
    % MODELO 1         X1  w   ph  x2
    par_var = [0.02 0.01 0.02 0.01]; % parameter variation in percent related to nominal
    % MODELO 2  X1  ph  x2
    %par_var = [0.02 0.02 0.02]; % parameter variation in percent related to nominal
    
    xr = xnom.*(1+(par_var/100).*(rand(1,length(xnom))-0.5));
    
    %regenerates the signal with the uncertainties
        i = 1;
        %MODELO 1
        %Xm = xr(1); Wf(i) = xr(2); Ph(i) = xr(3); KaS(i) = xr(4);
        %MODELO 2
        Xm = xr(1); Wf(i) = Wf; Ph(i) = xr(2); KaS(i) = xr(3);
        
        Ain = zeros(length(Xm),length(t));
        % Amplitude Step: applied after time passes 0
        Ain(i,:) = Xm(i);
        Ain(i,t >= 0) = Ain(i,t >= 0) * (1 + KxS(i));
        %Phase step
        Theta(i,:) = (Wf(i)*t) ...                         % Fundamental
                         + Ph(i);               % phase shift
        Theta(i,t >= 0) = Theta(i,t >= 0) + (KaS(i) * pi/180);
        cSignal = (Ain.*exp(-1i.*Theta));
        Signal = real(cSignal) + Aruido*(rand(1,length(t))-0.5);
        
        err = @(x) (Signal - f(x));
    y0 = f(x0);
    % figure
    % plot(t,Signal - y0,'r')
    % title('Point by point estimation error: Signal - y')
    
    % Levenberg-marquardt from optimization toolbox
    tol = 1e-7;
    OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','OptimalityTolerance',tol);
    OPTIONS.StepTolerance = 1e-12;
    OPTIONS.Display = 'off';
    [X,RESNORM,RESIDUAL,exitflag,output,lambda,jacobian] = lsqnonlin(err,x0,[],[],OPTIONS);
    Y = f(X);
    ferrHz(n) = X(2)/(2*pi) - F1;
    errors(n,:) = (xr - X)./xr;
    n
end

ERR_MAX = max(abs(errors))*100'   %erros maximos em %



figure
for k = 1:2
    subplot(2,2,k); hold on;
    histogram(errors(:,k));
end

MEAN_ERR = mean(errors)
STDEV_ERR = std(errors)
%hold on; plot(X(2),X(1),'rx', x0(2), x0(1), 'ro')

%freq_err = ERR_MAX(2)*F1/(2*pi)
ferrHz_max = max(abs(ferrHz))
    