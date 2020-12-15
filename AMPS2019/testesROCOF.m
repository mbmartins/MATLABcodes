%testes para estimação de fase, frequência e ROCOF
clear all; close all; clc
SNR = 90;
%angulo phi_0, fixo ou tabela
Pss = 0;
%fixed parameters
F0 = 60.0;
F1 = F0;
tau1 = 0.5;  % in [%] of the time window
tau2 = 2;
KaS = 10.0; %[degrees]
KxS = -0.1; % [relative step]
KRS = 0; %[Hz/s] %size of frequency step
Fs = 4800;
NCycles = 6;
NSamples = floor(NCycles*Fs/F0);

br = 0.05; brn = br*NSamples;
br_mask = floor((br)*NSamples):floor((1-br)*NSamples);

Signal = SigGEN2(F0,F1,Fs,Pss,NCycles,tau1,tau2,SNR,KaS, KxS,KRS);
z=hilbert(Signal');
theta_i = unwrap(angle(z));
df=gradient(theta_i);% Hilbert estimate of the instantaneous frequency of z
fest = median(df(brn:end-brn)); %frequency estimate
fest_hz = fest*Fs/(2*pi);
ROCOF_i = gradient(df);
ROCOF_est = median(ROCOF_i);

plot(Signal);
figure
plot(df(br_mask)*Fs/(2*pi)); title('estimated Fi [Hz]')
figure; plot(theta_i); title('estimated instant phase')
figure;plot((ROCOF_i(br_mask))); title('estimated ROCOF_i')

