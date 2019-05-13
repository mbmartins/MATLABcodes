%generation of performance comparison
clear all; close all; clc

%SNR = 30;
SNR =   [60 55  50  45  40  35 30];
HE_eps =  zeros(1,length(SNR));
NLHE_eps =  zeros(1,length(SNR));


for s = 1:length(SNR)
    eps_NLHE = 0;
    eps_HE = 0;
    %fixed parameters
    f0 = 60;
    tau1 = 0.5;  % in [%] of the time window
    SAG_cycles = 10; %duration of SAG
    KaS = 0; %[degrees]
    KxS = -0.1; % [relative step]
    Ps = 90;
    %MC for HE
    MCruns = 10000;
    for r = 1:MCruns  %taking the eps errors for tau1 only
        [tau_error_HE(r,:)] = HE_estimator(SNR(s),KxS,KaS,Ps,tau1,SAG_cycles);
        [FE(r), tau_error_NLHE(r,:)] = NLHE_estimator(SNR(s),KxS,KaS,Ps,tau1,SAG_cycles);
        eps_HE = eps_HE + (tau_error_HE(r,1)>2 || isnan(tau_error_HE(r,1)));
        eps_NLHE = eps_NLHE + (tau_error_NLHE(r,1)>2 || isnan(tau_error_NLHE(r,1)));
        
    end
    FE_std(s) = std(FE);
    FE_mean(s) = mean(FE);
    NLHE_eps(s) = eps_NLHE*100/MCruns; % in [%]
    HE_eps(s) = eps_HE*100/MCruns; % in [%]
end




figure
plot(SNR,HE_eps,'ko-',SNR,NLHE_eps,'bx--');
title('Estimation of \tau')
xlabel('SNR [dB]'); 
ylabel('N_e/N_{MC}'); legend('HE','NLHE')

figure
semilogy(SNR,FE_std,'ko-');
title('Estimation of frequency - standard deviation of FE')
xlabel('SNR [dB]'); 
ylabel('\sigma_{FE} [Hz]')

figure
histfit(FE) %do último valor de SNR