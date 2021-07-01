clear all; close all; clc;

load('waveforms.mat')

Fs = 4800;
Pin = 0;
dt = 1/Fs;
analysis_cycles = 6;
F0 = 60; T = 1/60;
samples_cycle = T/dt;
n_window_samples = samples_cycle*analysis_cycles;
fbin = Fs/n_window_samples;
freqs = [0:fbin:Fs/2-fbin];

%amplitude do sinal (Vp)
Xm = [0.9899 0.8 0 0.4 9.9985*0.8/12 9.9985*0.8/12 9.9985*0.8/12 8*0.8/12 0]';

wave1 = wav1.AmplitudePlot0;
wave2 = wav2.AmplitudePlot0;
wave3 = wav3.AmplitudePlot0;
wave4 = wav4.AmplitudePlot0;
wave5 = wav5.AmplitudePlot0;
wave6 = wav6.AmplitudePlot0; %amplif x10, 120:0.8
wave7 = wav7.AmplitudePlot0;
wave8 = wav8.AmplitudePlot0;
wave9 = wav9.AmplitudePlot0;

NSamples = length(wave1);
Nwindows = NSamples/n_window_samples;

ini = 1; %initial sample
%sig1 = wave1(ini:2^12);
sig1 = wave1(ini:4800);
t_window = 0:dt:dt*n_window_samples;

[snr_sig1, noise_pow1]  = snr(sig1)
snr(sig1,Fs,6)

sig1_model = Xm(1)*cos(2*pi*F0*t_window+Pin);
sig1_em = VSemulator(sig1_model,16,snr_sig1);

figure
[snr_sig1_em, noise_pow1_em]  = snr(sig1_em)
snr(sig1_em,Fs,6)

fft_sig1 = fft(sig1);
Xfft_sig1 = abs(fft_sig1)';

figure
semilogy(freqs,Xfft_sig1(1:length(freqs)))




sig2 = wave2(1:n_window_samples);
snr_sig2 = snr(sig2)

sig3 = wave3(1:n_window_samples); %zero signal
snr_sig3 = snr(sig3)
var_sig3 = std(sig3)


sig4 = wave4(1:n_window_samples);
snr_sig4 = snr(sig4)

sig5 = wave5(1:n_window_samples);
snr_sig5 = snr(sig5)


sig6 = wave6(1:n_window_samples);
snr_sig6 = snr(sig6)


sig7 = wave7(1:n_window_samples);
snr_sig7 = snr(sig7)

sig8 = wave8(1:n_window_samples);
snr_sig8 = snr(sig8)

sig9 = wave9(1:n_window_samples); %sinal com nivel zero
var_sig9 = std(sig9)
snr_sig9 = snr(sig9)

var_sig9/var_sig3