% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

%signal generation
F0 = 60; %Hz nominal
F1 = 60.0; %Hz fundamental
SampleRate = 4800; %Hz
dt = 1/SampleRate;
AnalysisCycles = 6;
NSamples = floor(AnalysisCycles*SampleRate/F0);
n = -NSamples/2:(NSamples/2-1); %discrete time vector
tau_pp = 0.1; % relative time of step in percent of total time 
tau_0 = (tau_pp - 0.5)*NSamples; %discrete time displacement
n = n - tau_0;
t = n*dt; %time vector
Vm = 100; %70*sqrt(2);
Ps = 0; %phase in degrees

% Phase in radians
Ph = Ps*pi/180;

KaS = 0;   % IEEE Std phase (angle) step index: 10 degrees
KxS = 0.1;   % magnitude step index: 0.1 
Wf = 2*pi*F1;  % fundamental frequency

Xm = Vm; %for now, single phase; TODO: 6-channels

%%%% MONTE CARLO
%%%% Signal parameters estimation
% Phase step
% Model: f(x) = x1*cos(w*t + phi + x2*u(t - tau))
u = zeros(length(Xm),length(t));
tau = 0;
u(length(Xm),t >= tau) = u(length(Xm),t >= tau) + 1;
f = @(x) x(1)*cos(x(2)*t + x(3) + x(4)*(pi/180)*u);
err = @(x) (Signal - f(x)).^2;
%first estimates for x
SNR = 99.5; %dB SNR = 20 log_10 Asinal/Aruido => Aruido = Asinal/10^(SNR/20)
Aruido = Vm/10^(SNR/20);

F1 = 60;
xnom1 = [Vm 2*pi*F1 Ph KaS];
x = xnom1;

xnom = xnom1;
Niter = 1000
x0 = xnom;  %x0 is fixed - first guess are the nominal values
for n = 1:Niter;
    %first guess
    %          X1  w   ph  x2
    par_var = [0 0.01 0 0]; % parameter variation in percent related to nominal
    xr = xnom1.*(1+(par_var/100).*(rand(1,length(xnom1))-0.5));
    
    %regenerates the signal with the uncertainties
        i = 1;
        Xm = xr(1); Wf(i) = xr(2); Ph(i) = xr(3); KaS(i) = xr(4);
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
        err = @(x) (Signal - f(x)).^2;
    y0 = f(x0);
    %Hilbert filter
    br = floor(0.05*NSamples); % 5% of NSamples are taken off at the beggining and end
    z=hilbert(Signal');  % calculates the analytic signal associated with Signal
    Psi = unwrap(angle(z));
    f_hi=unwrap(angle(z(2:end,:).*conj(z(1:end-1,:))));  % Hilbert estimate of the instantaneous frequency of Signal
    freq=f_hi-median(f_hi(br:end-br));
    %Ain = (Ain - mean(Ain))./abs(Ain);
    [ifmax, imax] = max(abs(freq(br:end-br)));
    imax_th = imax + br - 1;

    %estimation of tau
    if imax_th>0
        tau_e = imax_th*dt
        tau = (tau_0 + NSamples/2)*dt
    else
        tau_e = 0;
        tau = 0;
        imax_th = size(t,2)+1;
    end

    tau_error = (tau_e - tau)

    %splittig Signal in two
    imax_r = imax + br - 1;
    Signal_1 = Signal(1:imax_r-1);
    Signal_2 = Signal(imax_r+2:end);  %excluding samples around tau
    t_1 = t(1:imax_r-1);
    t_2 = t(imax_r+2:end);
    plot(t_1,Signal_1,t_2,Signal_2)

    %LSE of parameters 
    SignalParams = [F1 0 0 0 0 0 0 0];
    DelayCorr = 0;
    MagCorr = 1;

    [Synx1,Freq1,ROCOF1] = SteadyStateFit(	SignalParams, ...
        DelayCorr, ...
        MagCorr, ...
        F0, ...
        AnalysisCycles, ...
        SampleRate, ...
        Signal_1);
    [Synx2,Freq2,ROCOF2] = SteadyStateFit(	SignalParams, ...
        DelayCorr, ...
        MagCorr, ...
        F0, ...
        AnalysisCycles, ...
        SampleRate, ...
        Signal_2);

    fr(n) = xr(2)/(2*pi);
    ferr1(n) = (fr(n) - Freq1)*100;  %in [%]
    ferr2(n) = (fr(n) - Freq2)*100;  %in [%]
       
    %show iteration
    n
end

ferr1_max = max(abs(ferr1))*100'
ferr2_max = max(abs(ferr2))*100'

%
% ERR_MAX = max(abs(errors))*100'   %erros maximos em %
% 
% figure
% for k = 1:2
%     subplot(2,2,k); hold on;
%     histogram(errors(:,k));
% end
% 
% MEAN_ERR = mean(errors)
% STDEV_ERR = std(errors)




