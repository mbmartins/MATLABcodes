%generation of performance comparison
clear all; close all; clc

SNR = 50;

%angulo phi_0, fixo ou tabela
Pss = 0; %pior caso para salto mag phi = 0

%fixed parameters
F0 = 60.0;
F1 = 60.0;
Fs = 4800;
NCycles = 6;
T = NCycles/F0;
tau = 0.5;  % in [%] of the time window
tau2 = 10;
SAG_cycles = 10; %duration of SAG 
nbits = 16;

% %CASE 1 - salto de mag
h_x = 0.0; % [relative step]
h_a = 8.0; %[degrees]
%h_f = -1.0; %[Hz]
h_f = -0.; %[Hz]

%limiar para detector com PATV

% 
% parametro para PATV_HE
%lambda_f = 0.5:0.5:7;
lambda_f = 2.5;
lambda_m = 2;

Lf = 4.0e-3;
Lm = 1.9e-3;

%lambda = 0.11;
for k = 1:length(lambda_f); %para d=0

%lambda = 1.; %para d=1;

%maximo erro de tau toleravel em dt
max_dt = 2;
max_dt_7 = 7;
max_dt_8 = 8;

% em phi = 0, temos uma maior distribuição dos erros entao max_dt tem que
% ser em torno de 8
% essa dispersão diminui a medida em que phi aumenta

MCruns = 10000;

        %MC engine
        for r = 1:MCruns  %taking the eps errors for tau1 only
            % Fixed phi_0, from table
            Ps = Pss;
            
            Signal = SigGEN2(F0,F1,Fs,Ps,NCycles,tau,tau2,SNR,h_a,h_x,h_f); 
            [tau_ep,dmax_p] = HD_PATV_estimator_optim(Signal,lambda_m,lambda_f(k),Lm,Lf); 

            %tau_e é dado em dt
            %tau_e = tau_e/Fs; dt=1/Fs;
            NSamples = length(Signal);
            tau_error_HD_PATV(r) = (NSamples*tau_ep/T - floor(tau*NSamples)); %erro desta rodada de MC em dt
            
        
        %performance para PATV
         crit_tau_HD_PATV(r) = abs(tau_error_HD_PATV(r))> max_dt;
         crit_tau_HD_PATV_7(r) = abs(tau_error_HD_PATV(r))> max_dt_7;
         crit_tau_HD_PATV_8(r) = abs(tau_error_HD_PATV(r))> max_dt_8;
         crit_nan_HD_PATV(r) = isnan(tau_error_HD_PATV(r));
         crit_HD_PATV(r) = (crit_tau_HD_PATV(r) || crit_nan_HD_PATV(r));
         crit_HD_PATV_7(r) = (crit_tau_HD_PATV_7(r) || crit_nan_HD_PATV(r));
         crit_HD_PATV_8(r) = (crit_tau_HD_PATV_8(r) || crit_nan_HD_PATV(r));
        
        end
    false_pos_HD_PATV(k) = sum(crit_tau_HD_PATV)*100/MCruns;
    false_pos_HD_PATV_7(k) = sum(crit_tau_HD_PATV_7)*100/MCruns;
    false_pos_HD_PATV_8(k) = sum(crit_tau_HD_PATV_8)*100/MCruns;
    false_neg_HD_PATV(k) = sum(crit_nan_HD_PATV)*100/MCruns;

    std_tau_e_HD_PATV(k) = std(tau_error_HD_PATV);

    eps_HD_PATV(k) = sum(crit_HD_PATV)*100/MCruns; % in [%]; 
    eps_HD_PATV_7(k) = sum(crit_HD_PATV_7)*100/MCruns; % in [%]; 
    eps_HD_PATV_8(k) = sum(crit_HD_PATV_8)*100/MCruns; % in [%]; 
        
    message = "lambda:" + lambda_f(k)
    dmax_p_min(k) = min(dmax_p);
    dmax_p_mean(k) = mean(dmax_p);
end

%save('otimiza_lambda_f.mat')    

%limiar_p_mean = mean(limiar_p)

% figure(2)
% plot(lambda_f,eps_HD_PATV);


%graficos e analises antigas do HD e HD-PATV, verificar se pode ser util
%dmax_min = min(dmax)
% dmax_mean = mean(dmax)
% dmax_std = std(dmax)
%createfigure1(SNR, HE_eps')
%createfigure2(SNR, HE_eps')