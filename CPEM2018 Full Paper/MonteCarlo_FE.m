%monte carlo FE
clear all; close all; clc

SNRset = 50:10:100
MCruns = 500;
KxS = 0.1;
KaS = 0;
Ps = 0;
F1 = 60.0; %Hz

FE_LM_hist = zeros(length(SNRset),MCruns);
FE_SS_hist = zeros(length(SNRset),MCruns);
TVE_LM_hist = zeros(length(SNRset),MCruns);
TVE_SS_hist = zeros(length(SNRset),MCruns);

for snr = 1:length(SNRset)
    SNR = SNRset(snr);
    
    for k = 1:MCruns
        [FE_LM, FE_SS, TVE_LM, TVE_SS] = LM_comparison_FE(SNR,F1,Ps, KxS,KaS);
        FE_LM_hist(snr,k) = FE_LM;
        FE_SS_hist(snr,k) = FE_SS;
        TVE_LM_hist(snr,k) = TVE_LM;
        TVE_SS_hist(snr,k) = TVE_SS;
    end

% subplot(2,1,1)
    % hist(FE_LM_hist)
    % title('FE LM distribution')
    % subplot(2,1,2)
    % hist(FE_SS_hist)
    % title('FE SS distribution')

    FE_LM_mean(snr) = mean(FE_LM_hist(snr,:));
    FE_LM_std(snr) = std(FE_LM_hist(snr,:));
    FE_SS_mean(snr) = mean(FE_SS_hist(snr,:));
    FE_SS_std(snr) = std(FE_SS_hist(snr,:));

    TVE_LM_std(snr) = std(TVE_LM_hist(snr,:));
    TVE_LM_mean(snr) = mean(TVE_LM_hist(snr,:));
    TVE_SS_std(snr) = std(TVE_SS_hist(snr,:));
    TVE_SS_mean(snr) = mean(TVE_SS_hist(snr,:));

end

figure
semilogy(SNRset,FE_LM_std,'ro:',SNRset,FE_SS_std,'kx--')
ylabel('FE [Hz]'); xlabel('SNR dB');
legend('\sigma_{LM}','\sigma_{SS}')
title('FE accuracy: standard deviations, @F_1 = 60.0 Hz')

figure
semilogy(SNRset,TVE_LM_std,'ro:',SNRset,TVE_SS_std,'kx--')
ylabel('\sigma TVE [%]'); xlabel('SNR dB');
legend('\sigma TVE_{LM}','\sigma TVE_{DSS}')
title('TVE accuracy: standard deviations, @F_1 = 60.0 Hz')
%ylim([2.06e-3 2.12e-3])
%xlim([59.9 60.1])

figure
semilogy(SNRset,TVE_LM_mean,'ro:',SNRset,TVE_SS_mean,'kx--')
ylabel('TVE [%]'); xlabel('SNR dB');
legend('TVE_{LM}','TVE_{SS}')
title('TVE accuracy: average, @F_1 = 60.0 Hz')
% ylim([2.06e-3 2.12e-3])
% xlim([59.9 60.1])

% 
figure
semilogy(SNRset,abs(FE_LM_mean),'ro:',SNRset,abs(FE_SS_mean),'kx--')
ylabel('FE [Hz]'); xlabel('SNR dB');
legend('abs mean_{LM}','abs mean_{SS}')
title('FE accuracy: average, @F_1 = 60.0 Hz')

figure
subplot(2,1,1)
hist(FE_LM_hist(2,:),10);title('FE_{LM} @ SNR = 60dB')
xlim([-6e-4 6e-4])
subplot(2,1,2)
hist(FE_SS_hist(2,:),10); title('FE_{SS} @ SNR = 60dB')
xlim([-6e-4 6e-4])

figure
subplot(2,1,1)
hist(TVE_LM_hist(2,:));title('TVE_{LM} @SNR = 60dB')
%xlim([-8e-6 8e-6])
subplot(2,1,2)
hist(TVE_SS_hist(2,:)); title('TVE_{SS} @SNR = 60dB')
%xlim([-8e-6 8e-6])