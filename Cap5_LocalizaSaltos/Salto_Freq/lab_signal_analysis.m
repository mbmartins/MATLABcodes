% Analise das medições feitas para o Cap6
clear all; close all; clc;

load("meds_salto_freq.mat");

% ler ondas
% caso hf = 1 Hz
% wav_1 = KfSm1KaS10KmSm01Fs4800.AmplitudePlot0;
% wav_2 = KfSm1KaS10KmSm01Fs4800.AmplitudePlot1;
% wav_3 = KfSm1KaS10KmSm01Fs4800.AmplitudePlot2;
%ini_sample = 215; end_sample = ini_sample + 480 -1;

% caso hf = -0,5 Hz - não mudou???
% wav_1 = KfSm05KaS10KmSm01Fs4800.AmplitudePlot0;
% wav_2 = KfSm05KaS10KmSm01Fs4800.AmplitudePlot1;
% wav_3 = KfSm05KaS10KmSm01Fs4800.AmplitudePlot2;
% ini_sample = 1535-241; end_sample = ini_sample + 480 -1;

% caso hf = -1 Hz - não mudou???
wav_1 = KfSm1KaS10KmSm01Fs4800.AmplitudePlot0;
wav_2 = KfSm1KaS10KmSm01Fs4800.AmplitudePlot1;
wav_3 = KfSm1KaS10KmSm01Fs4800.AmplitudePlot2;
ini_sample = 456-241; end_sample = ini_sample + 480 -1;
plot(wav_1);hold on;
plot(wav_2);
plot(wav_3);

%wav 1 é a que contem o salto de frequencia

% caso hf = -1
wav_1 = newkm13.AmplitudePlot0;
wav_2 = newkm13.AmplitudePlot1;
wav_3 = newkm13.AmplitudePlot2;
ini_sample = 907-241; end_sample = ini_sample + 480 -1;
plot(wav_1);hold on;
plot(wav_2);
plot(wav_3);

% pegar janela de interesse
wav1 = wav_1(ini_sample:end_sample);
wav2 = wav_2(ini_sample:end_sample);
wav3 = wav_3(ini_sample:end_sample);
figure
plot(wav1);hold on;
plot(wav2);
plot(wav3);

%medicao de tau
kr = 3;
[tau_FD,a,b] = FD_estimator(wav1,kr);
lambda = 0.14; Lr = 1e-7;
[tau_FD_PATV,a,b] = FD_PATV_estimator(wav1,Lr,lambda);

tau_FD_error = tau_FD - 241
tau_FD_PATV_error = tau_FD_PATV - 241

