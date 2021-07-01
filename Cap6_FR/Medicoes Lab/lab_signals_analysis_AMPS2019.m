% Analise das medições feitas no artigo AMPS2019

clear all; close all; clc;

Fs = 4800; dt = 1/Fs; F0 = 60;
% th_gmi = 1e-3; th_fi = 1e-4; %thresholds for detection
% br = 0.02; %percentage of samples to be ignored in detection

wav_HN = load("meds_07_06.mat");
%Vp = 0.8Vp; sag duration 1 cycles; mag step -0.1
wav_1 = wav_HN.MAGONLYPS90;
wav_2 = wav_HN.PHONLYPS0;
wav_3 = wav_HN.BOTHPS0;

n1 = wav_1.TimePlot0;
n2 = wav_2.TimePlot0;
n3 = wav_3.TimePlot0;

ini_sample = 59*480 + 3*80 + 1; end_sample = ini_sample + 480;
wav1 = wav_1.AmplitudePlot0(n1 >= ini_sample & n1 < end_sample);
wav2 = wav_2.AmplitudePlot0(n2 >= ini_sample & n2 < end_sample);
wav3 = wav_3.AmplitudePlot0(n3 >= ini_sample & n3 < end_sample);

%plot(wav3)
SNR1 = snr(wav_1.AmplitudePlot0(1:4800))
SNR2 = snr(wav_2.AmplitudePlot0(1:4800))
SNR3 = snr(wav_3.AmplitudePlot0(1:4800))

%TODO:
% - adaptar o tau_estimator para estimação sem o PATV
% - ajustar os thresholds e parametros iguais aos da simulação
% - tentar medir os 3 casos simulados

F0 = 60;
km1 = 3; kf1 = 10;
km2 = 10; kf2 = 3;
km3 = 3; kf3 = 3;
[tau1,fest1] = HE_tau(km1,kf1,wav1,F0,Fs);
[tau2,fest2] = HE_tau(km2,kf2,wav2,F0,Fs);
[tau3,fest3] = HE_tau(km3,kf3,wav3,F0,Fs);

tau1n = tau1/dt
tau2n = tau2/dt
tau3n = tau3/dt

%[tau1,tau2,fest] = tau_estimator(wav1,F0,Fs)

FE1 = (fest1 - F0)/F0
FE2 = (fest2 - F0)/F0
FE3 = (fest3 - F0)/F0

% ---- HD-PATV
lambda_a = 2;
lambda_theta = 2.5;
[tau1_P,fest1_P, phi_01] = HD_PATV_estimator(wav1,lambda_a,lambda_theta);
[tau2_P,fest2_P, phi_02] = HD_PATV_estimator(wav2,lambda_a,lambda_theta);
[tau3_P,fest3_P, phi_03] = HD_PATV_estimator(wav3,lambda_a,lambda_theta);

tau1nP = tau1_P/dt
tau2nP = tau2_P/dt
tau3nP = tau3_P/dt

%[tau1,tau2,fest] = tau_estimator(wav1,F0,Fs)

FE1P = (fest1_P - F0)/F0
FE2P = (fest2_P - F0)/F0
FE3P = (fest3_P - F0)/F0

%consultar nas tabelas
corr(1) = 0.0;  %caso 1
corr(2) = 2.983345E-03; %caso 2, para tau ~ 0.5, phi_0 = 4 graus ~ 0;
corr(3) =  2.983345E-03;;  % caso 3 do AMPS...

FE1P_c = FE1P - corr(1)
FE2P_c = FE2P - corr(2)
FE3P_c = FE3P - corr(3)

% fazer medição de hf em wav1 e wav2
[wav1_f1,wav1_f2,wav1_f_r,wav1_h_f] = MedSF(wav1, tau1n, Fs)
[wav1_f1p,wav1_f2p,wav1_f_rp,wav1_h_fp] = MedSF_PATV(wav1, tau1n, Fs)

[wav2_f1,wav2_f2,wav2_f_r,wav2_h_f] = MedSF(wav2, tau2n, Fs)
[wav2_f1p,wav2_f2p,wav2_f_rp,wav2_h_fp] = MedSF_PATV(wav2, tau2n, Fs)

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


