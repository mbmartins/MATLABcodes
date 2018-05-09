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
SNR = 90; %dB SNR = 20 log_10 Asinal/Aruido => Aruido = Asinal/10^(SNR/20)
Aruido = Vm/10^(SNR/20);
Signal = real(cSignal) + Aruido*(rand(1,length(t))-0.5);

%we should filter the signal...
% N_FIR = 100;        % FIR filter order
% Fp  = 500;       % 500 Hz passband-edge frequency
% Fs  = SampleRate;       % 96 kHz sampling frequency
% Rp  = 0.00057565; % Corresponds to 0.01 dB peak-to-peak ripple
% Rst = 1e-4;       % Corresponds to 80 dB stopband attenuation
% 
% eqnum = firceqrip(N_FIR,Fp/(Fs/2),[Rp Rst],'passedge'); % eqnum = vec of coeffs
% lowpassFIR = dsp.FIRFilter('Numerator',eqnum); %or eqNum200 or numMinOrder
% %fvtool(lowpassFIR,'Fs',Fs,'Color','White')
% Sig_filt = lowpassFIR(Signal);
% figure; plot(Sig_filt,'r');
% legend('Filtered Signal');

%%%% Estimation of tau
%addpath('.\common')
%load filter_PB20_20Hz_5th % b20 a20 (elliptic lowpass filter with 20Hz passband and zero at 60 Hz)

%Hilbert filter
br = floor(0.05*NSamples); % 5% of NSamples are taken off at the beggining and end
z=hilbert(Signal');  % calculates the analytic signal associated with Signal
Psi = unwrap(angle(z));

f_hi=unwrap(angle(z(2:end,:).*conj(z(1:end-1,:))));  % Hilbert estimate of the instantaneous frequency of Signal
f=f_hi-median(f_hi(br:end-br));

%Ain = (Ain - mean(Ain))./abs(Ain);

%Total Variation of Psi
% m = length(Psi);
% D = zeros(m,m); 
% for i = 1:m-1
%    D(i,i) = -1;
%    D(i,i+1) = 1;
% end
% D(m,m) = 1;
% 
% TV = sum(abs(D*Psi));
% figure; 
% plot(t,TV,t,Psi); title('Total variation of \Psi')
% legend('TV(\Psi)','\Psi')


%threshold
th = 0.01;

[ifmax, imax] = max(abs(f(br:end-br)));
imax_th = find(abs(f(br:end-br))>th,1);
imax_th = imax_th + br - 1;

if imax_th>0
    tau_e = imax_th*dt
    tau = (tau_0 + NSamples/2)*dt
else
    tau_e = 0;
    tau = 0;
    imax_th = size(t,2)+1;
end

tau_error = (tau_e - tau)

%trying to demodulate signal with nominal values
% for phase
demod = 1; %cos(KaS*pi/180);
if imax_th>0
    Sig_demod = [Signal(1:imax_th) Signal(imax_th+1:end)+demod];
else 
    Sig_demod = Signal;
end
plot(t,Signal)
hold on; plot(t,Sig_demod); legend('Signal','Demodulated Signal')
%now, recalculate hilbert with demodulated signal
z=hilbert(Sig_demod');  % calculates the analytic signal associated with Signal
Psi = unwrap(angle(z));
f_hi=unwrap(angle(z(2:end,:).*conj(z(1:end-1,:))));  % Hilbert estimate of the instantaneous frequency of Signal
f=f_hi-median(f_hi(br:end-br));


%splittig Psi in two
Psi_1 = Psi(br:imax_th-1-br);
t_1 = t(br:imax_th-1-br);
Psi_2 = Psi(imax_th+br:end-br);
t_2 = t(imax_th+br:end-br);

%estimation of angular frequency (w = 2*pi*f) using hilbert
% w is the slope of the linear curve

%1st degree model matrix
%calculating 2freqs
X1 = [t_1' ones(size(t_1))']; y = Psi_1;
beta = (X1'*X1)\X1'*y;
fest_1 = beta(1)/(2*pi)
ferr_1 = 100*(fest_1 - F1)/F1
res_1 = (X1*beta - y);

X2 = [t_2' ones(size(t_2))']; y = Psi_2;
beta = (X2'*X2)\X2'*y;
fest_2 = beta(1)/(2*pi)
ferr_2 = 100*(fest_2 - F1)/F1
res_2 = (X2*beta - y);

fest = mean([fest_1 fest_2])
ferr = 100*(fest - F1)/F1

figure
subplot(3,1,1); plot(t_1,Psi_1,'.',t_2,Psi_2,'.'); title('\Psi')
subplot(3,1,2); plot(t_1,res_1,'.',t_2,res_2,'.'); title('residuals = (X*\beta - y)')
subplot(3,1,3); plot(t(2:end),abs(f_hi)); title('Instantaneous frequency')

%Other approach
%calculating 1freq
Psi = [Psi_1; Psi_2]; y = Psi;
t_m = [t_1 t_2]';
X = [t_m'; ones(size(t_m))']';
beta = (X'*X)\X'*y;
fest_m = beta(1)/(2*pi)
ferr_m = 100*(fest_m - F1)/F1
res_m = (X*beta - y);
figure
subplot(3,1,1); plot(t_m,Psi); title('\Psi')
subplot(3,1,2); plot(t_m,res_m,'.'); title('residuals = (X*\beta - y)')
subplot(3,1,3); plot(t(2:end),abs(f_hi)); title('Instantaneous frequency')


