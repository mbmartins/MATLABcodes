% Estimacao de frequencia e rocof 
%
% analise geral, com phi0 e tau aleatorios

clear all; close all; clc
SNR = 60;
%fixed parameters
F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency
k_a = 0.0; %[degrees]
k_x = 0.0; % [relative magnitude step]
k_f = 1.; %[Hz] %size of frequency step
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
% [FEraw_zero,fE1raw_zero,fE2raw_zero,kfEraw_zero] = MC_estimation_geral(MCiter,F0,F1,Fs,NCycles,SNR,k_a, k_x,k_f,nbits);
% 
% %salto fase
% k_a = 10; k_x = 0; k_f = 0;
% [FEraw_fase,fE1raw_fase,fE2raw_fase,kfEraw_fase] = MC_estimation_geral(MCiter,F0,F1,Fs,NCycles,SNR,k_a, k_x,k_f,nbits);
% %salto mag
% k_a = 0; k_x = 0.1; k_f = 0;
% [FEraw_mag,fE1raw_mag,fE2raw_mag,kfEraw_mag] = MC_estimation_geral(MCiter,F0,F1,Fs,NCycles,SNR,k_a, k_x,k_f,nbits);
%salto freq
k_a = 10; k_x = 0; k_f = 1;
[FEraw_freq,fE1raw_freq,fE2raw_freq,kfEraw_freq] = MC_estimation_geral(MCiter,F0,F1,Fs,NCycles,SNR,k_a, k_x,k_f,nbits);

save('boxplot_geral_media_Ftau60')
run('analise_geral_media_mediana')



% estimadores original:
% EF1,2,3: mean nos 3
% EF4: mean nos 3, mas parece ser melhor c mediana
% EF5: media de f1 e f2, fr por tau, funciona melhor c a media
% EF6: media de f1 e f2, fr por tau, funciona melhor c a mediana