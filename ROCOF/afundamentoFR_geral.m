%testes para estimação de fase, frequência e ROCOF
clear all; close all; clc
SNR = 60;
%fixed parameters
F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency
k_a = 10.0; %[degrees]
k_x = 0.0; % [relative magnitude step]
k_f = 0.; %[Hz] %size of frequency step
Fs =4800;
NCycles = 6;
T = NCycles/F0;
NSamples = floor(NCycles*Fs/F0);
nbits = 16;

MCiter = 5000;
phi_n = 360*rand(1,MCiter); %distribuicao de phi_0
tau_vec = 0.1 + (0.8)*rand(1,MCiter); %distribuicao de tau

% afundamento com duração de 0 a 0.8T, sem variacao de frequencia
k_a = 10; k_x = -0.1; k_f = 0.0;
duracao = 0.1;%*rand(1,MCiter);
tau_vec(2,:) = tau_vec(1,:) + duracao; %rand(1,MCiter); %distribuição de tau2;
%estimando frequencia com estimador para 1 salto somente
[FEraw_sag,fE1raw_sag,fE2raw_sag,kfEraw_sag, riraw_sag, draw_sag, tau_nraw_sag] = MC_estimation_ROCOF(MCiter,F0,F1,Fs,phi_n,NCycles,tau_vec,SNR,k_a, k_x,k_f,nbits);
kfest_sag = kfEraw_sag + k_f;

% afundamento com duração de 0.1T, sem variacao de frequencia
k_a = 30; k_x = -0.3; k_f = 0.0;
%estimando frequencia com estimador para 1 salto somente
[FEraw_sag30,fE1raw_sag30,fE2raw_sag30,kfEraw_sag30, riraw_sag30, draw_sag30, tau_nraw_sag30] = MC_estimation_ROCOF(MCiter,F0,F1,Fs,phi_n,NCycles,tau_vec,SNR,k_a, k_x,k_f,nbits);
kfest_sag30 = kfEraw_sag30 + k_f;


save('afundamento_geral')
run('analise_geral_afundamento')