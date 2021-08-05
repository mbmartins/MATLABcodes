%generation of performance comparison
clear all; close all; clc

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
h_f = -1; %[Hz]

%limiar para detector com PATV
Lr= 1e-7;
% 

% parametro para PATV_HE
lambda = 0.02:0.001:0.2;
for k = 1:length(lambda); %para d=0

%lambda = 1.; %para d=1;

%maximo erro de tau toleravel em dt
max_dt = 2;
% em phi = 0, temos uma maior distribuição dos erros entao max_dt tem que
% ser em torno de 8
% essa dispersão diminui a medida em que phi aumenta

MCruns = 1000;

        %MC engine
        for r = 1:MCruns  %taking the eps errors for tau1 only
            % Fixed phi_0, from table
            Ps = Pss;
            
            Signal = SigGEN2(F0,F1,Fs,Ps,NCycles,tau,tau2,SNR,h_a,h_x,h_f,nbits);
             
            [tau_ep,dmax_p(r),limiar_p(r)] = FD_PATV_estimator(Signal,Lr,lambda(k)); 

            %tau_e é dado em dt
            %tau_e = tau_e/Fs; dt=1/Fs;
            NSamples = length(Signal);
            tau_error_FD_PATV(r) = (tau_ep - floor(tau*NSamples)); %erro desta rodada de MC
            
        
        %performance para PATV
         crit_tau_FD_PATV(r) = abs(tau_error_FD_PATV(r))> max_dt;
         crit_nan_FD_PATV(r) = isnan(tau_error_FD_PATV(r));
         crit_FD_PATV(r) = (crit_tau_FD_PATV(r) || crit_nan_FD_PATV(r));
        
        end
        false_pos_FD_PATV(k) = sum(crit_tau_FD_PATV)*100/MCruns;
        false_neg_FD_PATV(k) = sum(crit_nan_FD_PATV)*100/MCruns;
        
        std_tau_e_FD_PATV(k) = std(tau_error_FD_PATV);
        
        eps_FD_PATV(k) = sum(crit_FD_PATV)*100/MCruns; % in [%]; 
    message = "lambda:" + lambda(k)
    dmax_p_min(k) = min(dmax_p);
    end

%limiar_p_mean = mean(limiar_p)

semilogy(lambda,dmax_p_min);


beep

%graficos e analises antigas do HD e HD-PATV, verificar se pode ser util
%dmax_min = min(dmax)
% dmax_mean = mean(dmax)
% dmax_std = std(dmax)
%createfigure1(SNR, HE_eps')
%createfigure2(SNR, HE_eps')