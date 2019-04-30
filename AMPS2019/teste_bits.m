%teste de precisão

clear all; close all; clc;

t = 0.01:0.01:10;
f = 1;
sinal = sin(2*pi*f*t) + 0.0000*randn(1,length(t));

digitsOld = digits(6);

sig = vpa(sinal)

sinal_d = double(sig)
figure
plot(t,sinal_d)
figure
snr(sinal_d)