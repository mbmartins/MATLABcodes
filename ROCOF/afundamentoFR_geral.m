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

% MC loop
%baseline
% k_a = 0; k_x = 0; k_f = 0;
% [FEraw_zero,fE1raw_zero,fE2raw_zero,kfEraw_zero] = MC_estimation_geral_afund(MCiter,F0,F1,Fs,NCycles,SNR,k_a, k_x,k_f,nbits);

%salto fase + salto mag combinados, frequencia constante
k_a = 10; k_x = 0.1; k_f = 0;
[FEraw_fase,fE1raw_fase,fE2raw_fase,kfEraw_fase] = MC_estimation_geral_afund(MCiter,F0,F1,Fs,NCycles,SNR,k_a, k_x,k_f,nbits);

save('afundamento_geral')
run('analise_geral_afundamento')