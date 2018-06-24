% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

%Monte carlo for tau_error
for k = 1:1000
%signal generation
F0 = 60; %Hz
F1 = 60; %Hz
UF1 = 0.04; %uncertainty of frequency in %
F1 = F1 + F1*2*UF1/100*(rand-0.5);

SampleRate = 5000; %Hz
dt = 1/SampleRate;
AnalysisCycles = 6;
NSamples = floor(AnalysisCycles*SampleRate/F0);
n = -NSamples/2:(NSamples/2-1); %discrete time vector
tau_pp_set = 0.1:0.1:0.9; % relative time of step in percent of total time 
tau_pp = tau_pp_set(randi(9));
tau_0 = (tau_pp - 0.5)*NSamples; %discrete time displacement
n = n - tau_0;
t = n*dt; %time vector
Vm = 0.9; %70*sqrt(2);
UVm = 1;
Vm = Vm + Vm*2*UVm/100*(rand-0.5);

Ps = 180; %phase in degrees
UPs = 100;
Ps = Ps + Ps*2*UPs/100*(rand-0.5);

Ps_deg(k) = Ps;

% Phase in radians
Ph = Ps*pi/180;

KaS = 0;   % IEEE Std phase (angle) step index: 10 degrees
UKa = 1;
KaS = KaS + KaS*2*UKa/100*(rand-0.5);

KxS = 0.1;   % magnitude step index: 0.1 
UKx = 1;
KxS = KxS + KxS*2*UKx/100*(rand-0.5);

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
SNR = 65; %dB SNR = 20 log_10 Asinal/Aruido => Aruido = Asinal/10^(SNR/20)
Aruido = Vm/10^(SNR/20);

  
    Signal = real(cSignal) + Aruido*(rand(1,length(t))-0.5);
    %plot(Signal)
    %%%% Estimation of tau
    br = 0.05*NSamples; % 5% of NSamples are taken off at the beggining and end
    z=hilbert(Signal');  % calculates the analytic signal associated with Signal
    f=angle(z(2:end,:).*conj(z(1:end-1,:)));  % Hilbert estimate of the instantaneous frequency of xds_f
    f=f-median(f(br:end-br));
    Ain = (Ain - mean(Ain))./abs(Ain);
    nn = 1:(NSamples-1);
    %figure
    %plot(nn,abs(f)); title('Instantaneous frequency (abs)')
    % figure
    % plot(f.^2); title('IF ^2');

    [ifmax, imax] = max(abs(f(br:end-br)));

    tau_e = (br + imax - 1)*dt;
    tau = (tau_0 + NSamples/2)*dt;

    tau_error(k) = tau_e - tau;
end

max_dt = max(tau_error)/dt
min_dt = min(tau_error)/dt

% figure
% subplot(2,1,1)
% plot(nn,abs(f)); title('Instantaneous frequency (abs)')
% subplot(2,1,2)
% plot(t,Signal); title('Signal with Magnitude Step')