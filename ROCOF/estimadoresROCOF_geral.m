%testes para estimação de fase, frequência e ROCOF

clear all; close all; clc
SNR = 60;
%fixed parameters
F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency
Fs =4800;
NCycles = 6;
T = NCycles/F0;
NSamples = floor(NCycles*Fs/F0);
nbits = 16;

MCiter = 2;
%phi_n = 360*rand(1,MCiter); %distribuicao de phi_0

%phi_0 fixo em zero:
phi_n = 0*ones(1,MCiter);
%tau_vec = 0.1 + (0.8)*rand(1,MCiter); %distribuicao de tau

%tau fixo em 0.5T
tau_vec = 0.5*ones(1,MCiter);;

% MC loop
for k = 1:MCiter
%baseline
k_a = 0; k_x = 0; k_f = 0;
 [FEraw_zero,fE1raw_zero,fE2raw_zero,kfEraw_zero, riraw_zero, draw_zero] = MC_estimation_ROCOF(MCiter,F0,F1,Fs,phi_n(k),NCycles,tau_vec(k),SNR,k_a, k_x,k_f,nbits);
%salto fase
k_a = 10; k_x = 0; k_f = 0;
[FEraw_fase,fE1raw_fase,fE2raw_fase,kfEraw_fase, riraw_fase, draw_fase] = MC_estimation_ROCOF(MCiter,F0,F1,Fs,phi_n(k),NCycles,tau_vec(k),SNR,k_a, k_x,k_f,nbits);
%salto mag
k_a = 0; k_x = 0.1; k_f = 0;
[FEraw_mag,fE1raw_mag,fE2raw_mag,kfEraw_mag, riraw_mag, draw_mag] = MC_estimation_ROCOF(MCiter,F0,F1,Fs,phi_n(k),NCycles,tau_vec(k),SNR,k_a, k_x,k_f,nbits);
%salto freq
k_a = 0; k_x = 0; k_f = 1;
[FEraw_freq,fE1raw_freq,fE2raw_freq,kfEraw_freq, riraw_freq, draw_freq] = MC_estimation_ROCOF(MCiter,F0,F1,Fs,phi_n(k),NCycles,tau_vec(k),SNR,k_a, k_x,k_f,nbits);
end
save('boxplot_ROCOF_geral')
run('analise_ROCOF_geral')