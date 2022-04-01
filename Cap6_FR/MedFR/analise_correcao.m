close all; clc; clear all;



%resultados anteriores
% freq_corr(1) = 0.00298276140145513; %phi0 = 0
% freq_corr(2) = 0.00296409387030716; %phi0 = 15
% freq_corr(3) = 0.00296598890104022; %phi0 = 30
% freq_corr(4) = 0.00306010061155422; %phi0 = 45
% freq_corr(5) = 0.00323290947446844; %phi0 = 60
% freq_corr(6) = 0.00338467357633335; %phi0 = 75
% freq_corr(7) = 0.00342652604283979; %phi0 = 90
% 
% phi_0 = [0:15:90];

% gradiente
% dfdphi = diff(freq_corr)/15;

% Passo 1 - fazer função para correção sistemática

N = 1000;
lambda = 2.5; F0 = 60; F1 = 60; Fs = 4800; NCycles = 6; tau1 = 0.5; tau2 = 1.;
SNR = 60; ha = 10; hm = 0; hf = 0; nbits = 16;
dP = 5;
Ps = -5:dP:180;
for k = 1:length(Ps)
    % MC loop
    for nmc = 1:N
        Signal = SigGEN2(F0,F1,Fs,Ps(k),NCycles,tau1,tau2,SNR, ha, hm,hf,nbits);
        [f_r_est(nmc), phi_0_est(nmc)] = MedFR_PATV(Signal,lambda);
        FE = (f_r_est(k) - F1)/F1;
    end
    Cf(k) = mean(FE);
    phi_0_mean(k) = mean(phi_0_est);
    phi_0std(k) = std(phi_0_est);
end

% Passo 2 - avaliar robustez da estimação de Cf em função de phi_0
% dCf/dphi
Ephi = Ps - phi_0_mean
dCf = gradient(Cf)/dP;

save('analise_correcao.mat') 

