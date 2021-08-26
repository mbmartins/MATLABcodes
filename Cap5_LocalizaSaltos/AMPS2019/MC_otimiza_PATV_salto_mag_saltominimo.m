%generation of performance comparison
clear all; close all; clc

% roteiro: determinar primeiro o minimo salto que se deseja detectar, com
% criterios fisicos: hf = 0,5Hz. Daí, calcular o Lf desse salto minimo em funçao de lambda.
% Escolher depois um lambda que otimize o desempenho para o salto nominal.

%OBS: resultados pessimos para o FD, não vale a pena insistir muito
% apostar no PATV

SNR = 60;

%angulo phi_0, fixo ou tabela
%Pss = 0;  %pior caso para salto de freq é phi0 = 0;
Pss = 45;  %pior caso para determinar Lr parece ser 45;
%Pss = 0:15:90;

%fixed parameters
F0 = 60.0;
F1 = 60.0;
Fs = 4800;
NCycles = 6;
tau = 0.5;  % in [%] of the time window
tau2 = 10;
SAG_cycles = 10; %duration of SAG 
nbits = 16;

% %CASE 3 - salto de frequencia
h_a = 0.0; %[degrees]
h_x = -0.; % [relative step]
%h_f = -1.0; %[Hz]
h_f = [-0.5 -0.7 -0.9]; %[Hz]

%limiar para detector com PATV
Lr= 0;
% 

% parametro para PATV_HE
%lambda = 0.02:0.01:0.4;
lambda = 0.11;

for k = 1:length(lambda); %para d=0

%lambda = 1.; %para d=1;

%maximo erro de tau toleravel em dt
max_dt = 2;
max_dt_7 = 7;
max_dt_8 = 8;

% em phi = 0, temos uma maior distribuição dos erros entao max_dt tem que
% ser em torno de 8
% essa dispersão diminui a medida em que phi aumenta

MCruns = 1000;

        %MC engine
        for r = 1:MCruns  %taking the eps errors for tau1 only
            % Fixed phi_0, from table
            Ps = Pss;
            
            Signal1 = SigGEN2(F0,F1,Fs,Ps,NCycles,tau,tau2,SNR,h_a,h_x,h_f(1),nbits);
            [tau_ep,dmax_p1(r),limiar_p1(r)] = FD_PATV_estimator(Signal1,Lr,lambda(k)); 

            Signal2 = SigGEN2(F0,F1,Fs,Ps,NCycles,tau,tau2,SNR,h_a,h_x,h_f(2),nbits);
            [tau_ep,dmax_p2(r),limiar_p2(r)] = FD_PATV_estimator(Signal2,Lr,lambda(k)); 
            
            Signal3 = SigGEN2(F0,F1,Fs,Ps,NCycles,tau,tau2,SNR,h_a,h_x,h_f(3),nbits);
            [tau_ep,dmax_p3(r),limiar_p3(r)] = FD_PATV_estimator(Signal3,Lr,lambda(k)); 
            
            
            %tau_e é dado em dt
            %tau_e = tau_e/Fs; dt=1/Fs;
            NSamples = length(Signal1);
            tau_error_FD_PATV(r) = (tau_ep - floor(tau*NSamples)); %erro desta rodada de MC
            
        
        %performance para PATV
         crit_tau_FD_PATV(r) = abs(tau_error_FD_PATV(r))> max_dt;
         crit_tau_FD_PATV_7(r) = abs(tau_error_FD_PATV(r))> max_dt_7;
         crit_tau_FD_PATV_8(r) = abs(tau_error_FD_PATV(r))> max_dt_8;
         crit_nan_FD_PATV(r) = isnan(tau_error_FD_PATV(r));
         crit_FD_PATV(r) = (crit_tau_FD_PATV(r) || crit_nan_FD_PATV(r));
         crit_FD_PATV_7(r) = (crit_tau_FD_PATV_7(r) || crit_nan_FD_PATV(r));
         crit_FD_PATV_8(r) = (crit_tau_FD_PATV_8(r) || crit_nan_FD_PATV(r));
        
        end
        false_pos_FD_PATV(k) = sum(crit_tau_FD_PATV)*100/MCruns;
        false_pos_FD_PATV_7(k) = sum(crit_tau_FD_PATV_7)*100/MCruns;
        false_pos_FD_PATV_8(k) = sum(crit_tau_FD_PATV_8)*100/MCruns;
        false_neg_FD_PATV(k) = sum(crit_nan_FD_PATV)*100/MCruns;
        
        std_tau_e_FD_PATV(k) = std(tau_error_FD_PATV);
        
        eps_FD_PATV(k) = sum(crit_FD_PATV)*100/MCruns; % in [%]; 
        eps_FD_PATV_7(k) = sum(crit_FD_PATV_7)*100/MCruns; % in [%]; 
        eps_FD_PATV_8(k) = sum(crit_FD_PATV_8)*100/MCruns; % in [%]; 
        
    message = "lambda:" + lambda(k)
    dmax_p_min1(k) = min(dmax_p1);
    dmax_p_mean1(k) = mean(dmax_p1);
    dmax_p_min2(k) = min(dmax_p2);
    dmax_p_mean2(k) = mean(dmax_p2);
    dmax_p_min3(k) = min(dmax_p3);
    dmax_p_mean3(k) = mean(dmax_p3);
end

%save('otimiza_lambda_r.mat')    

%limiar_p_mean = mean(limiar_p)

% figure(2)
% semilogy(lambda,dmax_p_min);


%graficos e analises antigas do HD e HD-PATV, verificar se pode ser util
%dmax_min = min(dmax)
% dmax_mean = mean(dmax)
% dmax_std = std(dmax)
%createfigure1(SNR, HE_eps')
%createfigure2(SNR, HE_eps')