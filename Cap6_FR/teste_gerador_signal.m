%testes gerador
close all; clear all; clc;
F0 = 60;
F1 = 58;
Fs = 4800;
phi_n = 0;
NCycles = 6;
tau1 = 0.5;
tau2 = 10*tau1;
SNR = 600;
KaS = 0;
KxS = 0;
KfS = 1;
nbits = 16;

Signal1 = SigGEN2(F0,F1,Fs,phi_n,NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits);
F1 = 60;
Signal2 = SigGEN2(F0,F1,Fs,phi_n,NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits);

%subplot(1,2,1)
plot(Signal1); grid on; hold on;
%subplot(122)
plot(Signal2); grid on;
