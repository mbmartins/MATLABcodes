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

MCiter = 5000;
phi_n = 360*rand(1,MCiter); %distribuicao de phi_0
%phi_n = 0*ones(1,MCiter);
tau_vec(1,:) = 0.1 + (0.8)*rand(1,MCiter); %distribuicao de tau1
tau_vec(2,:) = 1; %rand(1,MCiter); %distribuição de tau2;
%tau_vec = 0.5*ones(1,MCiter);;

% ------ MC loop
%baseline
k_a = 0; k_x = 0; k_f = 0;
 [FEraw_zero,fE1raw_zero,fE2raw_zero,kfEraw_zero, riraw_zero, draw_zero, tau_nraw_zero] = MC_estimation_ROCOF(MCiter,F0,F1,Fs,phi_n,NCycles,tau_vec,SNR,k_a, k_x,k_f,nbits);
kfest_zero = kfEraw_zero + k_f;

%salto fase
k_a = 10; k_x = 0; k_f = 0;
[FEraw_fase,fE1raw_fase,fE2raw_fase,kfEraw_fase, riraw_fase, draw_fase, tau_nraw_fase] = MC_estimation_ROCOF(MCiter,F0,F1,Fs,phi_n,NCycles,tau_vec,SNR,k_a, k_x,k_f,nbits);
kfest_fase = kfEraw_fase + k_f;

%salto mag
k_a = 0; k_x = 0.1; k_f = 0;
[FEraw_mag,fE1raw_mag,fE2raw_mag,kfEraw_mag, riraw_mag, draw_mag, tau_nraw_mag] = MC_estimation_ROCOF(MCiter,F0,F1,Fs,phi_n,NCycles,tau_vec,SNR,k_a, k_x,k_f,nbits);
kfest_mag = kfEraw_mag + k_f;

%salto freq
k_a = 0; k_x = 0; k_f = 1.0;
[FEraw_freq,fE1raw_freq,fE2raw_freq,kfEraw_freq, riraw_freq, draw_freq, tau_nraw_freq] = MC_estimation_ROCOF(MCiter,F0,F1,Fs,phi_n,NCycles,tau_vec,SNR,k_a, k_x,k_f,nbits);
kfest_freq = kfEraw_freq + k_f;

% salto fasemag
%salto freq
k_a = 10; k_x = 0.1; k_f = 0.0;
[FEraw_fm,fE1raw_fm,fE2raw_fm,kfEraw_fm, riraw_fm, draw_fm, tau_nraw_fm] = MC_estimation_ROCOF(MCiter,F0,F1,Fs,phi_n,NCycles,tau_vec,SNR,k_a, k_x,k_f,nbits);
kfest_fm = kfEraw_fm + k_f;

% afundamento com duração de 0.3T, sem variacao de frequencia
k_a = 10; k_x = 0.1; k_f = 0.0;
tau_vec(2,:) = tau_vec(1,:) + 0.3; %rand(1,MCiter); %distribuição de tau2;
[FEraw_sag,fE1raw_sag,fE2raw_sag,kfEraw_sag, riraw_sag, draw_sag, tau_nraw_sag] = MC_estimation_ROCOF(MCiter,F0,F1,Fs,phi_n,NCycles,tau_vec,SNR,k_a, k_x,k_f,nbits);
kfest_sag = kfEraw_sag + k_f;


% % ------ END MC LOOP 
% save('boxplot_ROCOF_geral_dmax.mat')
% run('analise_ROCOF_geral_dmax')

save('boxplot_ROCOF_geral_kfest.mat')
run('analise_ROCOF_geral_kf')