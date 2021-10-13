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
% wav_1 = KfSm1KaS10KmSm01Fs4800.AmplitudePlot0;
% wav_2 = KfSm1KaS10KmSm01Fs4800.AmplitudePlot1;
% wav_3 = KfSm1KaS10KmSm01Fs4800.AmplitudePlot2;
% ini_sample = 456-241; end_sample = ini_sample + 480 -1;
% plot(wav_1);hold on;
% plot(wav_2);
% plot(wav_3);

% caso hf = -1 Hz - não mudou???
% wav_1 = KfSm1KaS10KmSm01Fs4800.AmplitudePlot0;
% wav_2 = KfSm1KaS10KmSm01Fs4800.AmplitudePlot1;
% wav_3 = KfSm1KaS10KmSm01Fs4800.AmplitudePlot2;
% ini_sample = 456-241; end_sample = ini_sample + 480 -1;
% plot(wav_1);hold on;
% plot(wav_2);
% plot(wav_3);

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

%estimar tau
tau1n = 240;

% fazer medição de hf
[wav1_f1,wav1_f2,wav1_f_r,wav1_h_f] = MedSF(wav1, tau1n, Fs)
[wav1_f1p,wav1_f2p,wav1_f_rp,wav1_h_fp] = MedSF_PATV(wav1, tau1n, Fs)

% erros de fr
hf = -1;
Fref = (F0 + (hf + F0))/2;

FE_MedSF_caso3 = (wav1_f_r - Fref)
FE_MedSF_caso3 = (wav1_f_rp - Fref)

% erros de hf
hfE_MedSF_caso3 = wav1_h_f - hf
hfE_MedSF_caso3 = wav1_h_fp - hf



