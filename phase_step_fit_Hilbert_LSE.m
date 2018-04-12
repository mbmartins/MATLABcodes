% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

%signal generation
F0 = 60; %Hz nominal
F1 = 60; %Hz fundamental
SampleRate = 4800; %Hz
dt = 1/SampleRate;
AnalysisCycles = 6;
NSamples = floor(AnalysisCycles*SampleRate/F0);
n = -NSamples/2:(NSamples/2-1); %discrete time vector
tau_pp = 0.5; % relative time of step in percent of total time 
tau_0 = (tau_pp - 0.5)*NSamples; %discrete time displacement
n = n - tau_0;
t = n*dt; %time vector
Vm = 100; %70*sqrt(2);
Ps = 0; %phase in degrees

% Phase in radians
Ph = Ps*pi/180;

KaS = 10;   % IEEE Std phase (angle) step index: 10 degrees
KxS = 0;   % magnitude step index: 0.1 
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
SNR = 97; %dB SNR = 20 log_10 Asinal/Aruido => Aruido = Asinal/10^(SNR/20)
Aruido = 0*Vm/10^(SNR/20);
Signal = real(cSignal) + Aruido*(rand(1,length(t))-0.5);

%Hilbert filter
br = floor(0.05*NSamples); % 5% of NSamples are taken off at the beggining and end
z=hilbert(Signal');  % calculates the analytic signal associated with Signal
Psi = unwrap(angle(z));
f_hi=unwrap(angle(z(2:end,:).*conj(z(1:end-1,:))));  % Hilbert estimate of the instantaneous frequency of Signal
f=f_hi-median(f_hi(br:end-br));

%Ain = (Ain - mean(Ain))./abs(Ain);

[ifmax, imax] = max(abs(f(br:end-br)));
imax_th = imax;

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
imax_r = imax + br;
Signal_1 = Signal(1:imax_r);
Signal_2 = Signal(imax_r+1:end);
t_1 = t(1:imax_r);
t_2 = t(imax_r+1:end);
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


