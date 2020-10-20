%testes para estimação de fase, frequência e ROCOF

%variacoes interessantes:
% f1 = 59, 61, etc
% dois saltos
% encontrar um lambda otimo para o PATV

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
[FEraw,fE1raw,fE2raw,kfEraw] = MC_estimation_geral(MCiter,F0,F1,Fs,NCycles,SNR,k_a, k_x,k_f,nbits);
            
boxplot(FEraw')
xlabel('EF')
ylabel('FE')

save('salto_fase_geral')
% se desejar figuras rodar
%run('analise_lambda.m')