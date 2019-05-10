%generation of performance comparison
clear all; close all; clc

SNR = 60;
%SNR =   [60 55  50  45  40  35 30];
HE_eps =    [0  6.6 8.2 9   10  97 99];
NLHE_eps =  zeros(1,length(SNR));


for s = 1:length(SNR)
    eps = 0;
    %fixed parameters
    f0 = 60;
    tau1 = 0.5;  % in [%] of the time window
    SAG_cycles = 1; %duration of SAG
    KaS = 0; %[degrees]
    KxS = -0.1; % [relative step]
    Ps = 90;
    %MC for HE
    MCruns = 100;
    for r = 1:MCruns
        [FE(r), tau_error(r,:)] = NLHE_estimator(SNR(s),KxS,KaS,Ps,tau1,SAG_cycles);
        eps = eps + (tau_error(r,1)>3);
    end
    FE_std(s) = std(FE);
    FE_mean(s) = mean(FE);
    NLHE_eps(s) = eps*100/MCruns; % in [%]
end




figure
plot(SNR,HE_eps,'ko-',SNR,NLHE_eps,'bx-');
title('Estimation of \tau')
xlabel('SNR [dB]'); 
ylabel('N_e/N_{MC}')

figure
semilogy(SNR,FE_std,'ko-');
title('Estimation of frequency - standard deviation of FE')
xlabel('SNR [dB]'); 
ylabel('\sigma_{FE} [Hz]')

figure
histfit(FE) %do último valor de SNR