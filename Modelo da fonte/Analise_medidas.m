clear all;close all; clc;

medidas = load('medidas.mat')

%dados
ts = 7.81e-5;
fs = 1/ts;
amostras_ciclo = 128;

%janela de observação
ciclos = 14;
ciclo_inicio = 1;
ini = ciclo_inicio*amostras_ciclo+1; fin = ini + ciclos*amostras_ciclo - 1;

%leitura do sinal
sig1 = medidas.seno1.AmplitudeSLAVE(ini:fin);
%retirar nivel DC
%sig1 = sig1 - mean(sig1);
n1 = medidas.seno1.TimeSLAVE(ini:fin);
t1 = n1*ts;
N1 = size(n1,1)

%sinal emulado
f0 = 50.008; Vp = 5; SNR = 60; nbits = 12;
Input = 1*sin(2*pi*f0*t1)';
vout = Vp*VSemulator(Input,nbits,SNR);

% sig2 = medidas.seno1.AmplitudeMASTER;
% n2 = medidas.seno1.TimeMASTER;
% t2 = n2*ts;
% plot(t1,sig1,t2,sig2);

sig1_fft = fft(sig1);
sig2_fft = fft(vout);

%window_han = hanning(N1);
%sig1_w = sig1.*window_han;
%sig1_fftw = fft(sig1_w);
mag_sig1 = abs(sig1_fft(1:N1/2));
mag_sig2 = abs(sig2_fft(1:N1/2));
%mag_sig1w = abs(sig1_fftw(1:N1/2));
fbin = fs/N1;
freqs = ((1:N1/2)-1)*fbin;
semilogy(freqs,mag_sig1,freqs,mag_sig2)
legend('sinal amostrado','sinal modelo')


snr_amostrado = snr(sig1)
snr_emulado = snr(vout)

% figure
% subplot(2,1,1)
% snr(sig1, fs, 3)
% subplot(2,1,2)
% snr(vout, fs, 3)

% figure
% spectrogram(sig1)

