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
% - adaptar o tau_estimator para estima��o sem o PATV
% - ajustar os thresholds e parametros iguais aos da simula��o
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