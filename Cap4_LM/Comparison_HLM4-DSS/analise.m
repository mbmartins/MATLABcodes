clear all; close all; clc;

load("comparisonHLM4_DSS_SM10000.mat")

figure
subplot(121)
semilogy(SNRset,FE_LM_std,'ro-',SNRset,FE_SS_std,'kx-')
hold on; 
grid on;
% semilogy(SNRset,abs(FE_LM_mean),'ro--',SNRset,abs(FE_SS_mean),'kx--')
ylabel('FE [Hz]'); xlabel('SNR dB');
% lgd = legend('\sigma_{HLM4}','\sigma_{DSS}','\mu_{HLM4}','\mu_{DSS}')
lgd = legend('\sigma_{HLM4}','\sigma_{DSS}')
lgd.FontSize = 14;
%lgd.Location = 'Southoutside';
lgd.Orientation = 'Horizontal';
%title('FE - Salto de Magnitude')

subplot(122)
semilogy(SNRset,TVE_LM_std,'ro-',SNRset,TVE_SS_std,'kx-')
hold on; 
grid on;
% semilogy(SNRset,abs(FE_LM_mean),'ro--',SNRset,abs(FE_SS_mean),'kx--')
ylabel('TVE [%]'); xlabel('SNR dB');
% lgd = legend('\sigma_{HLM4}','\sigma_{DSS}','\mu_{HLM4}','\mu_{DSS}')
lgd = legend('\sigma_{HLM4}','\sigma_{DSS}')
lgd.FontSize = 14;
%lgd.Location = 'Southoutside';
lgd.Orientation = 'Horizontal';
%title('FE - Salto de Magnitude')

load("comparisonHLM4_DSS_SF1000.mat")

figure
subplot(121)
semilogy(SNRset,FE_LM_std,'ro-',SNRset,FE_SS_std,'kx-')
hold on; 
grid on;
% semilogy(SNRset,abs(FE_LM_mean),'ro--',SNRset,abs(FE_SS_mean),'kx--')
ylabel('FE [Hz]'); xlabel('SNR dB');
% lgd = legend('\sigma_{HLM4}','\sigma_{DSS}','\mu_{HLM4}','\mu_{DSS}')
lgd = legend('\sigma_{HLM4}','\sigma_{DSS}')
lgd.FontSize = 14;
%lgd.Location = 'Southoutside';
lgd.Orientation = 'Horizontal';
%title('FE - Salto de Magnitude')

subplot(122)
semilogy(SNRset,TVE_LM_std,'ro-',SNRset,TVE_SS_std,'kx-')
hold on; 
grid on;
% semilogy(SNRset,abs(FE_LM_mean),'ro--',SNRset,abs(FE_SS_mean),'kx--')
ylabel('TVE [%]'); xlabel('SNR dB');
% lgd = legend('\sigma_{HLM4}','\sigma_{DSS}','\mu_{HLM4}','\mu_{DSS}')
lgd = legend('\sigma_{HLM4}','\sigma_{DSS}')
lgd.FontSize = 14;
%lgd.Location = 'Southoutside';
lgd.Orientation = 'Horizontal';
%title('FE - Salto de Magnitude')
