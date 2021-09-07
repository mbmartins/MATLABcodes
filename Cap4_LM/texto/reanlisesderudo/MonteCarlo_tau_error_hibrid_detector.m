% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

UF1 = 0.05; %uncertainty of frequency in [%]
UVm = 1; %uncertainty of magnitude in [%]
UKa = 1; %uncertainty of Ka in [%]
UKx = 1; %uncertainty of Kx in [%]
UPs = 1; %uncertainty of phase in degrees
Pin = 90; %initial phase in degrees  %worst case is 90 degrees

SNR = 50;  %signal noise ratio in dB
rng('shuffle');

extremos = 0;

%Monte carlo for tau_error
for k = 1:1000

%signal generation
F0 = 60; %Hz
F1 = 60; %Hz
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
Vm = 0.9; 
Vm = Vm + Vm*2*UVm/100*(rand-0.5);

Ps = Pin; %phase in degrees
Ps = Ps + 2*UPs/100*(rand-0.5);

Ps_deg(k) = Ps;

% Phase in radians
Ph = Ps*pi/180;

KaS = 0;   % IEEE Std phase (angle) step index: 10 degrees
KaS = KaS + KaS*2*UKa/100*(rand-0.5);

KxS = 0.1;   % magnitude step index: 0.1 
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
 %dB SNR = 20 log_10 Asinal/Aruido => Aruido = Asinal/10^(SNR/20)
%Anoise = Vm/(sqrt(2))/10^(SNR/20);
rSignal = real(cSignal);
var_noise = ((std(rSignal))/(10^(SNR/20)))^2;
std_noise = (std(rSignal))/(10^(SNR/20));
noise = std_noise*randn(1,length(rSignal));  
%noise_un = std_noise*rand(1,length(rSignal));

    Signal =  rSignal + noise;
    %SNR_hist(k) = snr(Signal,noise);
    SNR_est(k)=10*log10(var(rSignal)/var(noise));
    %plot(Signal)
    %%%% Estimation of tau
    %br = 0.05*NSamples; % 5% of NSamples are taken off at the beggining and end
    br = 0.02*NSamples; % 5% of NSamples are taken off at the beggining and end
    z=hilbert(Signal');  % calculates the analytic signal associated with Signal
    %f=angle(z(2:end,:).*conj(z(1:end-1,:)));  % Hilbert estimate of the instantaneous frequency of xds_f
    f=gradient(phase(z));
    f=f-median(f(br:end-br));
    Ain = (Ain - mean(Ain))./abs(Ain);
    %nn = 1:(NSamples-1);
    nn = 1:(NSamples);
    
   
    [ifmax, imax] = max(abs(f(br:end-br)));
    limiar1=7*median(abs(f(br:end-br)));
    
    gmi=gradient(abs(z));
    gmi=gmi-median(gmi(br:end-br));
    [ifmax2, imax2] = max(abs(gmi(br:end-br)));
    limiar2=7*median(abs(gmi(br:end-br)));
    
    if ifmax(1)>limiar1,
        tau_e=(br + imax(1)-1)*dt;
    elseif ifmax2(1)>limiar2,
        tau_e=(br + imax2(1)-1)*dt;
    else
        tau_e=NaN;
    end

    %tau_e = (br + imax - 1)*dt;
    %tau_e = (br + imax(1)-1)*dt;
    tau = (tau_0 + NSamples/2)*dt;

    tau_error(k) = tau_e - tau;
    
    if abs(tau_error(k)/dt) > 2.1 
        extremos = extremos + 1;
        %figure
        %subplot(2,1,1)
        %plot(nn,abs(f));
        %xlabel('Samples', 'Fontsize',17);ylabel('$\mid{d[n]}\mid$','Interpreter','latex','Fontsize',17)
        %subplot(2,1,2)
        %plot(nn,Signal(2:500)); xlabel('Samples', 'Fontsize',17); ylabel('y(t) [V]', 'Fontsize',17)
    end 
end



max_dt = max(tau_error)/dt
min_dt = min(tau_error)/dt

figure
%hist((tau_error/dt),200)
hist((tau_error/dt),-250:250)
xlabel('\tau error [\Deltat]')
ylabel('Occurrences')

extremos
erro_tau=extremos*100/k
% figure
% subplot(2,1,1)
% plot(nn,abs(f)); title('Instantaneous frequency (abs)')
% subplot(2,1,2)
% plot(t,Signal); title('Signal with Magnitude Step')