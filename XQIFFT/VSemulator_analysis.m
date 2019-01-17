close all; clear all; clc;

%voltage source emulator
format long;

%signal
t = [0:1:5000];
fs = 5000;     
f1 = 60
Input1 = sin((2*pi*f1*t)/fs);

%parameters
bits = 12;   %33250 Agilent - 12 bits SNR -86dB
SNR = 60; %[dB]

quant=max(Input1)/(2^(bits-1)-1);
y=round(Input1/quant);
%signe=uint8((sign(y)'+1)/2)
%out=[signe dec2bin(abs(y),7)]  % The first bit represents the sign of the number

vout = y/max(abs(y));
figure 
title('Quantized signal')
plot(t,vout)

quantization_errors = Input1 - vout;
% figure
% hist(quantization_errors)
% title('Quantization errors histogram')

%spurious noise
    var_sig = std(Input1);
    eta = var_sig/10^(SNR/20); %eq (3) CPEM
    %vout2 = vout + eta*(randn(1,length(t))-0.5);
    vout2 = vout + normrnd(0,eta,[1,length(t)])

spurious_noise = vout2 - vout;
    
sqnr = 20*log10(2^bits)
sqnr2 = 6.02*bits
sqnr3 = 1.761 + sqnr2
snr_vout = snr(vout)
snr_vout2 = snr(vout2)

total_noise = vout2 - Input1;
figure 
subplot(3,1,1)
nbins = 20; dist = 'kernel';
histfit(quantization_errors,nbins,dist)
xlim([-3e-3 3e-3])
title('Quantization Noise')
subplot(3,1,2)
histfit(spurious_noise,nbins,'Normal')
title('Spurious Noise')
xlim([-3e-3 3e-3])
subplot(3,1,3)
nbins = 20; dist = 'Normal';
histfit(total_noise,nbins,dist)
title('Total Noise')
xlim([-3e-3 3e-3])



figure
x = [quantization_errors' spurious_noise' total_noise'];
boxplot(x,'Labels',{'Quantization','Spurious','Total'})
title('Comparison of noise components')

figure;
subplot(2,1,1)
normplot(spurious_noise); title('Spurious Noise - SNR = 60dB')
subplot(2,1,2)
normplot(total_noise); title('Total Noise')

var_quantization = std(quantization_errors)
var_noise = std(total_noise)
pd = fitdist(total_noise','Normal')

[hq,pq] = lillietest(quantization_errors)
% if hq == 0
%    fprint('Distribution passed goodness of fit test') 
% else
%    fprint('Distribution failed goodness of fit test') 
% end
[ht,pt] = lillietest(total_noise)

%incluir um passa baixa na saida 10MHz
% incluir ruido branco do nivel DC - gaussiano 60dB
% qual a influencia entao da quantizacao com esse ruido de 60dB???
