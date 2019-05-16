clear all; close all; clc;

Fs = 4800; dt = 1/Fs; F0 = 60;
% th_gmi = 1e-3; th_fi = 1e-4; %thresholds for detection
% br = 0.02; %percentage of samples to be ignored in detection

wav_HN = load("wav_HN.mat");
%Vp = 0.8Vp; sag duration 1 cycles; mag step -0.1
wav_aux = wav_HN.Vp081cycle;
n1 = wav_aux.TimePlot0;
ini_sample = 28600; end_sample = ini_sample + 480;
wav1 = wav_aux.AmplitudePlot0(n1 >= ini_sample & n1 < end_sample);

SNR = snr(wav_aux.AmplitudePlot0(1:2*4800))

[tau1,tau2,fest] = tau_estimator(wav1,F0,Fs)

FE = fest - F0
