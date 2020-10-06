%Investigação da estimação de fase incial
clear all; close all; clc
SNR = 90;
%angle phi_0 [degrees]
phi_0 = 0;
%fixed parameters
F0 = 10.0;
F1 = F0;
tau1 = 0.5;  % in proportion of the time window
tau2 = 2; % in proportion of the time window
KaS = 10.0; %[degrees]
KxS = 0.0; % [relative magnitude step]
KfS = 0.; %[Hz] %size of frequency step
Fs =4800;
NCycles = 6;
NSamples = floor(NCycles*Fs/F0);
tau_n = floor(tau1*NSamples);

Signal = SigGEN2(F0,F1,Fs,phi_0,NCycles,tau1,tau2,SNR,KaS, KxS,KfS);
z=hilbert(Signal');
Psi_i = unwrap(angle(z)); %[rad]
mz = abs(z);
%Psi_i = Psi_i.*mz./median(mz)
subplot(2,1,1)
plot(Psi_i*180/pi);
ylabel('Psi[n] [°]'); xlabel('Samples')
subplot(2,1,2)
plot(Psi_i(1:30)*180/pi);
ylabel('Psi[n] [°]'); xlabel('Samples')

df = gradient(Psi_i);
%df = diff(Psi_i);
df = df.*mz./median(mz);
% Psi_comp = cumsum(df);
% plot(Psi_i(1:10)); hold on; plot(Psi_comp(1:10))

phi_est = Psi_i(1)*180/pi
% phi_est_comp = Psi_comp(1)*180/pi