%generation of performance comparison
clear all; close all; clc

SNR = 60;
%SNR = 60:-10:30;
%SNR = 90:-10:30;

%angulo phi_0, fixo ou tabela
Pss = 0;
%Pss = 0:15:90;

%fixed parameters
f0 = 60.0;
tau1 = 0.5;  % in [%] of the time window
SAG_cycles = 10; %duration of SAG 

%CASE 1
%  KaS = 0.0; %[degrees]
%  KxS = -0.1; % [relative step]

% %CASE 2
KaS = 10.0; %[degrees]
KxS = -0.; % [relative step]

%CASE 3
%  KaS = 10.0; %[degrees]
%  KxS = -0.1; % [relative step]

%limiares para HE
km = 3;
kf = 3; % a partir de 8 praticamente fase não atua...

%limiares para PATV_HE
lambda_a = 2.;
%lambda_theta = 2.5;
lambda_theta = 2.5; %para freq, melhor lambda pequeno

MCruns = 1000;

for ps = 1:length(Pss) % loop for different initial phases
    for s = 1:length(SNR)  % loop for different SNRs
        %MC engine
        for r = 1:MCruns  %taking the eps errors for tau1 only
            % Fixed phi_0
            Ps = Pss(ps);
            
            % random phi_0 from uniform distribution
            % Ps_rand(r) = 90*rand(1,1);
            % Ps = Ps_rand(r);
            
            [tau_error_HE(r,:),FE_HE(r,:)] = HE_estimator(SNR(s),KxS,KaS,Ps,tau1,SAG_cycles,km,kf);
            [tau_error_PATV_HE(r,:),FE_PATV_HE(r,:),dmax(r,:)] = PATV_HE_estimator(SNR(s),KxS,KaS,Ps,tau1,SAG_cycles,lambda_a,lambda_theta);
            
        %Performance of tau estimation    
        crit_tau_HE(r) = abs(tau_error_HE(r,1))>2;
        crit_nan_HE(r) = isnan(tau_error_HE(r,1));
        crit_HE(r) = (crit_tau_HE(r) || crit_nan_HE(r));
        
         crit_tau_PATV_HE(r) = abs(tau_error_PATV_HE(r,1))>2;
         crit_nan_PATV_HE(r) = isnan(tau_error_PATV_HE(r,1));
         crit_PATV_HE(r) = (crit_tau_PATV_HE(r) || crit_nan_PATV_HE(r));
        
        end
        false_pos_HE(ps,s) = sum(crit_tau_HE)*100/MCruns;
        false_neg_HE(ps,s) = sum(crit_nan_HE)*100/MCruns;
         false_pos_PATV_HE(ps,s) = sum(crit_tau_PATV_HE)*100/MCruns;
         false_neg_PATV_HE(ps,s) = sum(crit_nan_PATV_HE)*100/MCruns;
        
        eps_HE(ps,s) = sum(crit_HE)*100/MCruns; % in [%];
        eps_PATV_HE(ps,s) = sum(crit_PATV_HE)*100/MCruns; % in [%];
        
        FE_std_HE(ps,s) = std(FE_HE)/f0;
        FE_mean_HE(ps,s) = mean(FE_HE)/f0;
        FE_box(s,:) = FE_HE/f0; %[Hz/Hz]

         FE_std_PATV_HE(ps,s) = std(FE_PATV_HE)/f0;
         FE_mean_PATV_HE(ps,s) = mean(FE_PATV_HE)/f0;
         FE_box(s,:) = FE_PATV_HE/f0; %[Hz/Hz]
        
    end
end


FE_std_HE
FE_std_PATV_HE
eps_HE
eps_PATV_HE
beep

dmax_min = min(dmax)

%FE_corr = FE_mean_PATV_HE - 0.002213
% dmax_mean = mean(dmax)
% dmax_std = std(dmax)
%createfigure1(SNR, HE_eps')
%createfigure2(SNR, HE_eps')

% boxplot(FE_box',SNR)
% xlabel('SNR [dB]')
% ylabel('FE [Hz/Hz]')
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