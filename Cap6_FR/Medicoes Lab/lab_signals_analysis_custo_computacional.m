% Analise das medições feitas no artigo AMPS2019

clear all; close all; clc;

Fs = 4800; dt = 1/Fs; F0 = 60;
% th_gmi = 1e-3; th_fi = 1e-4; %thresholds for detection
% br = 0.02; %percentage of samples to be ignored in detection
% 
% wav_HN = load("meds_07_06.mat");
% %Vp = 0.8Vp; sag duration 1 cycles; mag step -0.1
% wav_1 = wav_HN.MAGONLYPS90;
% wav_2 = wav_HN.PHONLYPS0;
% wav_3 = wav_HN.BOTHPS0;
% 
% n1 = wav_1.TimePlot0;
% n2 = wav_2.TimePlot0;
% n3 = wav_3.TimePlot0;
% 
% ini_sample = 59*480 + 3*80 + 1; end_sample = ini_sample + 480;
% wav1 = wav_1.AmplitudePlot0(n1 >= ini_sample & n1 < end_sample);
% wav2 = wav_2.AmplitudePlot0(n2 >= ini_sample & n2 < end_sample);
% wav3 = wav_3.AmplitudePlot0(n3 >= ini_sample & n3 < end_sample);
% 
% %plot(wav3)
% SNR1 = snr(wav_1.AmplitudePlot0(1:4800))
% SNR2 = snr(wav_2.AmplitudePlot0(1:4800))
% SNR3 = snr(wav_3.AmplitudePlot0(1:4800))

%TODO:
% - adaptar o tau_estimator para estimação sem o PATV
% - ajustar os thresholds e parametros iguais aos da simulação
% - tentar medir os 3 casos simulados

F0 = 60;
NMC = 10000;
N = 480;
tau = 0.1 + 0.8*rand(1,NMC);
tau_n = round(tau*N);
phi_0 = pi*rand(1,NMC);
n = 1:480;

hx = -0.1;
ha = 10;


%sinais
for nr=1:NMC
    u = n>tau_n(nr);
    %salto mag
    wav1(:,nr) = (1+ hx*u).*cos(2*pi*F0*n/Fs + phi_0(nr));
    %salto fase
    wav2(:,nr) = cos(2*pi*F0*n/Fs + phi_0(nr) + ha*u*pi/180);
end
% MedFR
tic
for nr = 1:NMC
[fest1] = MedFR(wav1(:,nr),Fs);
end
time_MedFR_wav1 = toc/NMC;
tic
for nr = 1:NMC
[fest2] = MedFR(wav2(:,nr),Fs);
end
time_MedFR_wav2 = toc/NMC;

% FE1 = (fest1 - F0)/F0
% FE2 = (fest2 - F0)/F0
% FE3 = (fest3 - F0)/F0

% ---- MedFR_PATV
lambda_a = 2;
lambda_theta = 2.5;
tic
for nr = 1:NMC
[fest1_P, phi_01] = MedFR_PATV(wav1(:,nr),lambda_theta,Fs);
end
time_MedFR_PATV_wav1 = toc/NMC;
tic
for nr = 1:NMC
[fest2_P, phi_02] = MedFR_PATV(wav2(:,nr),lambda_theta,Fs);
end 
time_MedFR_PATV_wav2 = toc/NMC;

%[tau1,tau2,fest] = tau_estimator(wav1,F0,Fs)

FE1P = (fest1_P - F0)/F0
FE2P = (fest2_P - F0)/F0

%consultar nas tabelas
corr(1) = 0.0;  %caso 1
corr(2) = 2.983345E-03; %caso 2, para tau ~ 0.5, phi_0 = 4 graus ~ 0;
corr(3) =  2.983345E-03;;  % caso 3 do AMPS...

FE1P_c = FE1P - corr(1)
FE2P_c = FE2P - corr(2)

% fazer medição de hf em wav1 e wav2
tau_n = 240; % o MedSF precisa da estimativa de taun
tic
for nr = 1:NMC
[wav1_f1,wav1_f2,wav1_f_r,wav1_h_f] = MedSF(wav1(:,nr), tau_n, Fs);
time_MedSF_wav1 = toc/NMC;
end
tic
for nr = 1:NMC
[wav1_f1p,wav1_f2p,wav1_f_rp,wav1_h_fp] = MedSF_PATV(wav1(:,nr), tau_n, Fs);
end
time_MedSF_PATV_wav1 = toc/NMC;

tic
for nr = 1:NMC
[wav2_f1,wav2_f2,wav2_f_r,wav2_h_f] = MedSF(wav2(:,nr), tau_n, Fs);
end
time_MedSF_wav2 = toc/NMC;

tic
for nr = 1:NMC
[wav2_f1p,wav2_f2p,wav2_f_rp,wav2_h_fp] = MedSF_PATV(wav2(:,nr), tau_n, Fs);
end
time_MedSF_PATV_wav2 = toc/NMC;

% erros de fr
FE_MedSF_caso1 = (wav1_f_r - F0)/F0
FE_MedSF_caso2 = (wav2_f_r - F0)/F0
FE_MedSF_caso1 = (wav1_f_rp - F0)/F0
FE_MedSF_caso2 = (wav2_f_rp - F0)/F0

% erros de hf
hfE_MedSF_caso1 = wav1_h_f
hfE_MedSF_caso2 = wav2_h_f
hfE_MedSF_caso1 = wav1_h_fp
hfE_MedSF_caso2 = wav2_h_fp


%medições de fr com LM
%salto mag: x0 = [Vm KxS 2*pi*F1 Ph];
F1 = F0; %pequeno erro de 3% na freq inicial
x1 = [1 -0.1 2*pi*F1 pi/2]; KaS = 0; KxS = -0.1;
tic
params = [0 0 0 0]; %alocacao de memoria
for nr = 1:NMC
[params] = LM_estimator(wav1(:,nr)',x1,tau_n, KaS, KxS, Fs);
end
time_LM_smag = toc/NMC;
fmed1 = params(3)/(2*pi)
ferror1 = (fmed1 - F0)/F0

% salto fase [Vm 2*pi*F1 Ph KaS];
x2 = [1 2*pi*F1 0 10]; KaS = 10; KxS = 0;
tic
params2 = [0 0 0 0];
for nr = 1:NMC
[params2] = LM_estimator(wav2(:,nr)',x2,tau_n, KaS, KxS, Fs);
end
time_LM_sfase = toc/NMC;
fmed_LM_sfase = params2(2)/(2*pi)
ferror_LM_sfase = (fmed_LM_sfase - F0)/F0

save('custo_computacional')
