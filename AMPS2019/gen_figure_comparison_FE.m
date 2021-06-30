%generation of performance comparison
clear all; close all; clc

%SNR = 60;
SNR = 60:-10:40;
%SNR = 90:-10:30;

%angulo phi_0, fixo ou tabela
%Pss = 0;
Pss = 0:30:90;

%fixed parameters
f0 = 60.;
tau1 = 0.5;  % in [%] of the time window
SAG_cycles = 10; %duration of SAG 

%CASE 1
  KaS = 0.0; %[degrees]
  KxS = -0.1; % [relative step]

% %CASE 2
%  KaS = 10.0; %[degrees]
%  KxS = -0.; % [relative step]

%CASE 3
%  KaS = 10.0; %[degrees]
%  KxS = -0.1; % [relative step]

%limiares para PATV_HE
lambda_a = 2.;
lambda_theta = 2.5;

%relatado no primeiro draft: lambda_a = .5; lambda_theta = 1.
lambda_a = 2.;
%lambda_theta = 2.5;
lambda_theta = 2.5; %para freq, melhor lambda pequeno

MCruns = 10000;

for ps = 1:length(Pss) % loop for different initial phases
    for s = 1:length(SNR)  % loop for different SNRs
        %MC engine
        for r = 1:MCruns  %taking the eps errors for tau1 only
            % Fixed phi_0
%             Ps = Pss(ps);
            
            % random phi_0 from uniform distribution
             Ps_rand(r) = 90*rand(1,1);
             Ps = Ps_rand(r);
            km = 3; kf=3;
            [tau_error_HE(r,:),FE_HE(r,:)] = HE_estimator(SNR(s),KxS,KaS,Ps,tau1,SAG_cycles,km,kf);
            [FE_PATV_HE(r,:),FE2_PATV_HE(r,:), FE3_PATV_HE(r,:)] = PATV_FE_estimator(SNR(s),KxS,KaS,Ps,tau1,SAG_cycles,lambda_a,lambda_theta);
        end
       
         FE_std_HE(ps,s) = std(FE_HE)/f0;
         FE_mean_HE(ps,s) = mean(FE_HE)/f0;
         FE_box(s,:) = FE_HE/f0; %[Hz/Hz]
         
         FE_std_PATV_HE(ps,s) = std(FE_PATV_HE)/f0;
         FE_mean_PATV_HE(ps,s) = mean(FE_PATV_HE)/f0;
         FE_box(s,:) = FE_PATV_HE/f0; %[Hz/Hz]
         
         FE2_std_PATV_HE(ps,s) = std(FE2_PATV_HE)/f0;
         FE2_mean_PATV_HE(ps,s) = mean(FE2_PATV_HE)/f0;
         FE2_box(s,:) = FE2_PATV_HE/f0; %[Hz/Hz]
         
         FE3_std_PATV_HE(ps,s) = std(FE3_PATV_HE)/f0;
         FE3_mean_PATV_HE(ps,s) = mean(FE3_PATV_HE)/f0;
         FE3_box(s,:) = FE3_PATV_HE/f0; %[Hz/Hz]
        
    end
end


beep

save('Case1_taufixo_psrand.mat')

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