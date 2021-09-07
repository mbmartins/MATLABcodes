% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox
% LM3 with fixed frequency
% Phase step
clear all; close all; clc;

% Initial parameters
F0 = 60; %Hz
F1 = 60; %Hz
SampleRate = 5000; %Hz
dt = 1/SampleRate;
AnalysisCycles = 6;
NSamples = floor(AnalysisCycles*SampleRate/F0);
n = -NSamples/2:(NSamples/2-1); %discrete time vector
SNR = 60; %dB SNR = 20 log_10 Asinal/Aruido => Aruido = Asinal/10^(SNR/20)
Vm = 1; %70*sqrt(2) =~ 100;
Xm = Vm;
Ps = -120; %phase in degrees
Ph = Ps*pi/180;% Phase in radians
KaS = -10;   % IEEE Std phase (angle) step index: 10 degrees
KxS = 0;   % magnitude step index: 0.1 
Wf = 2*pi*F1;  % fundamental frequency

 
for ti = 5
    % time vectors
    tau_pp = 0.1*ti; % relative time of step in percent of total time 
    tau_0 = (tau_pp - 0.5)*NSamples; %discrete time displacement
    n = -NSamples/2:(NSamples/2-1); %discrete time vector
    n = n - tau_0;
    t = n*dt; %time vector

    % step funtion 
    u = zeros(length(Xm),length(t));  tau = 0;
    u(length(Xm),t >= tau) = u(length(Xm),t >= tau) + 1;
    
    % Phase step Model: f(x) = x1*cos(w*t + phi + x2*u(t - tau)) 
    f = @(x) x(1)*cos(Wf*t + x(2) + 2*pi + x(3)*(pi/180)*u);
    xnom = [Vm Ph KaS];
    err = @(x) (Signal - f(x));

    % Nonlinear fit
    % Monte Carlo analysis
    Nruns = 1000;
    x0 = xnom;  %x0 is fixed - first guess are the nominal values
    for k = 1:Nruns
        %first guess
        k
        %uncertainties of parameters in signal generation (uniform
        %distribution)
        par_var = [1 1 1]; % parameter variation in percent related to nominal
        rng('shuffle'); rn = (rand(1,length(xnom))-0.5);
        xr = xnom.*(1+2*(par_var/100).*rn);
        Uwf = 0.000; % frequency uncertainty in [%] 
        Wf1 = Wf*(1 + (Uwf/100)*randn(1));
        %freq_rand = xr(3)/(2*pi)
        
        %uncertainties of tau estimation
        utau = 2;  %number of dts 
        u = zeros(length(Xm),length(t));
        tau = dt*randi([-utau utau],1);
        u(length(Xm),t >= tau) = u(length(Xm),t >= tau) + 1;
        
        %regenerates the signal with the uncertainties
        i = 1;
        Xm = xr(1); 
        Wf(i) = Wf1; 
        Ph(i) = xr(2);
        KxS_(i) = 0;
        KaS_(i) = xr(3);
        Ain = zeros(length(Xm),length(t));
        % Amplitude Step: applied after time passes 0
        Ain(i,:) = Xm(i);
        Ain(i,t >= 0) = Ain(i,t >= 0) * (1 + KxS_(i));
        %Phase step
        Theta(i,:) = (Wf(i)*t) ...                         % Fundamental
                         + Ph(i);               % phase shift
        Theta(i,t >= 0) = Theta(i,t >= 0) + (KaS_(i) * pi/180);
        cSignal = (Ain.*exp(-1i.*Theta));
        rSignal = real(cSignal);
        var_noise = ((std(rSignal))/(10^(SNR/20)))^2;
        std_noise = (std(rSignal))/(10^(SNR/20));
        noise = std_noise*randn(1,length(rSignal));
        Signal = rSignal + noise;
%        SNR_hist = snr(Signal,noise);

        err = @(x) (Signal - f(x));
        y0 = f(x0);
        x = xnom;

        % Levenberg-marquardt from optimization toolbox
        tol = 1e-7;
        OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','OptimalityTolerance',tol);
        OPTIONS.StepTolerance = 1e-12;
        OPTIONS.Display = 'iter-detailed';
        [X,RESNORM,RESIDUAL,exitflag,output,lambda, jacobian] = lsqnonlin(err,x0,[],[],OPTIONS);
        Y = f(X);
        
        %Intermediate phasor
        T = NSamples*dt;
        tau_est = tau_pp + tau/NSamples;
        Xe_r = xr(1);
        Xe = X(1);
        Phe_r = xr(2)*tau_pp + (xr(2) + xr(3)*(pi/180))*(1 - tau_pp);
        Phe = X(2)*tau_est + (X(2)+X(3)*(pi/180))*(1 - tau_est);
        
        Phasor_mag_error(k,:) = 100*(Xe - Xe_r)/Xe_r; % [%]
        Phasor_ph_error(k,:) = (Phe - Phe_r)*180/pi;   % [degrees]

    end

    Phasor_mag_errmean(ti,1) = mean(Phasor_mag_error);
    Phasor_ph_errmean(ti,1) = mean(Phasor_ph_error);
    Phasor_mag_errstdev(ti,1) = std(Phasor_mag_error);
    Phasor_ph_errstdev(ti,1) = std(Phasor_ph_error);        
       
end    

beep