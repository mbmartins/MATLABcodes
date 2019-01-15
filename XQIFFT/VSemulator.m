close all; clear all; clc;

%voltage source emulator
format long;
%voltages = [1.0 0. -0.99]
% 
% temp = voltages;
% mask = temp < 0;
% temp(mask) = 2^12 + temp(mask) ;
% binary = dec2bin(temp, 12);
% 
% temp_max = 2^12;
% vout = temp/temp_max

% outro exemplo
bits = 12;   %33250 Agilent - 12 bits SNR -86dB
t = [0:1:500];
fs = 5000;     
Input1 = sin((2*pi*50*t)/fs);
quant=max(Input1)/(2^(bits-1)-1);
y=round(Input1/quant);
%signe=uint8((sign(y)'+1)/2)
%out=[signe dec2bin(abs(y),7)]  % The first bit represents the sign of the number

vout = y/max(abs(y));
figure 
title('Quantized signal')
plot(t,vout)

quantization_errors = Input1 - vout;
figure
hist(quantization_errors)
title('Quantization errors histogram')

sqnr = 20*log10(2^bits)
sqnr2 = 6.02*bits
sqnr3 = 1.761 + sqnr2
snr_vout = snr(vout)

%incluir um passa baixa na saida 10MHz
% incluir ruido branco do nivel DC