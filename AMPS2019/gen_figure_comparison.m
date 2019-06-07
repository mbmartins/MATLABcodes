%generation of performance comparison
clear all; close all; clc

%SNR = 30;
%SNR = 60:-10:30;
SNR = 60:-10:30;
Pss = 0:15:90;

%angulo fixo
%Pss = 90;

%fixed parameters
f0 = 60;
tau1 = 0.5;  % in [%] of the time window
SAG_cycles = 10; %duration of SAG 
KaS = 10; %[degrees]
KxS = -0.0; % [relative step]

km = 10;
kf = 3; % a partir de 8 praticamente fase não atua...

HE_eps =  zeros(1,length(SNR));
NLHE_eps =  zeros(1,length(SNR));
MCruns = 10000;

for ps = 1:length(Pss)
    for s = 1:length(SNR)
        eps_NLHE = 0;
        eps_HE = 0;
        %MC for HE
        for r = 1:MCruns  %taking the eps errors for tau1 only
            Ps = Pss(ps);
            %Ps_rand(r) = 90*rand(1,1);
            %Ps = Ps_rand(r);
            [tau_error_HE(r,:),FE(r,:)] = HE_estimator(SNR(s),KxS,KaS,Ps,tau1,SAG_cycles,km,kf);
    %        [FE_NL(r), tau_error_NLHE(r,:)] = NLHE_estimator(SNR(s),KxS,KaS,Ps,tau1,SAG_cycles);
        crit_tau(r) = abs(tau_error_HE(r,1))>2;
        crit_nan(r) = isnan(tau_error_HE(r,1));
        crit_HE(r) = (crit_tau(r) || crit_nan(r));
    %        eps_NLHE = eps_NLHE + (tau_error_NLHE(r,1)>2 || isnan(tau_error_NLHE(r,1)));
        end
        false_pos(ps,s) = sum(crit_tau)*100/MCruns;
        false_neg(ps,s) = sum(crit_nan)*100/MCruns;
        eps_HE = sum(crit_HE);
        FE_std(ps,s) = std(FE)/f0;
        FE_mean(ps,s) = mean(FE)/f0;
        FE_box(s,:) = FE/f0; %[Hz/Hz]
        %NLHE_eps(s) = eps_NLHE*100/MCruns; % in [%]
        HE_eps(ps,s) = eps_HE*100/MCruns; % in [%]
    end
end

HE_eps
%createfigure1(SNR, HE_eps')
%createfigure2(SNR, HE_eps')

boxplot(FE_box',SNR)
xlabel('SNR [dB]')
ylabel('FE [Hz/Hz]')
% 
% figure
% plot(SNR,FE_mean, SNR,FE_mean + FE_std, SNR, FE_mean - FE_std)
% 
% 
% figure
% semilogy(SNR,FE_std,'ko-');
% % title('Estimation of frequency - standard deviation of FE')
% xlabel('SNR [dB]'); 
% ylabel('\sigma_{FE} [Hz]')
% 
%figure; histfit(FE_box(7,:))
% figure
% histfit(FE) %do último valor de SNR