%generation of performance comparison
clear all; close all; clc

%OBS: resultados pessimos, não vale a pena insistir muito
% apostar no PATV

SNR = 60;
%SNR = 60:-10:30;
%SNR = 90:-10:30;

%angulo phi_0, fixo ou tabela
%Pss = 45;  %pior caso para salto de freq é phi0 = 45;
Pss = 0:15:90;

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

%fator de multiplicação para limiares de detecção
kf = 10; % para salto de frequencia pq 10x??
%limiar para detector com PATV
Lf = 1e-13;
% 
% parametro para PATV_HE
lambda = 0.15; %para d=0
%lambda = 1.; %para d=1;

%maximo erro de tau toleravel em dt
max_dt = 8;
% em phi = 0, temos uma maior distribuição dos erros entao max_dt tem que
% ser em torno de 8
% essa dispersão diminui a medida em que phi aumenta

MCruns = 1000;

for ps = 1:length(Pss) % loop for different initial phases
    for s = 1:length(SNR)  % loop for different SNRs
        %MC engine
        for r = 1:MCruns  %taking the eps errors for tau1 only
            % Fixed phi_0, from table
            Ps = Pss(ps);
            
            % random phi_0 from uniform distribution
%              Ps_rand(r) = 90*rand(1,1);
%              Ps = Ps_rand(r);
            
             Signal = SigGEN2(F0,F1,Fs,Ps,NCycles,tau,tau2,SNR(s),h_a,h_x,h_f,nbits);
             
            [tau_e,dmax(r),limiar(r)] = FD_estimator(Signal,kf); %kf é o limiar de detecção
            [tau_ep,dmax_p(r),limiar_p(r)] = FD_PATV_estimator(Signal,Lf,lambda); %kf é o limiar de detecção

            %tau_e é dado em dt
            %tau_e = tau_e/Fs; dt=1/Fs;
            NSamples = length(Signal);
            tau_error_FD(r) = (tau_e - floor(tau*NSamples)); %erro desta rodada de MC
            tau_error_FD_PATV(r) = (tau_ep - floor(tau*NSamples)); %erro desta rodada de MC
            
        %Performance of tau estimation    
        crit_tau_FD(r) = abs(tau_error_FD(r))> max_dt;
        crit_nan_FD(r) = isnan(tau_error_FD(r));
        crit_FD(r) = (crit_tau_FD(r) || crit_nan_FD(r));
        
        %performance para PATV
         crit_tau_FD_PATV(r) = abs(tau_error_FD_PATV(r))> max_dt;
         crit_nan_FD_PATV(r) = isnan(tau_error_FD_PATV(r));
         crit_FD_PATV(r) = (crit_tau_FD_PATV(r) || crit_nan_FD_PATV(r));
        
        end
        false_pos_FD(ps,s) = sum(crit_tau_FD)*100/MCruns;
        false_neg_FD(ps,s) = sum(crit_nan_FD)*100/MCruns;
        false_pos_FD_PATV(ps,s) = sum(crit_tau_FD_PATV)*100/MCruns;
        false_neg_FD_PATV(ps,s) = sum(crit_nan_FD_PATV)*100/MCruns;
        
        eps_FD(ps,s) = sum(crit_FD)*100/MCruns % in [%];
        eps_FD_PATV(ps,s) = sum(crit_FD_PATV)*100/MCruns % in [%]; 
    message = "phi:"+Ps+"°,  SNR:"+SNR(s)+" dB, MCrun:"+r
    end
end

limiar_mean = mean(limiar)
dmax_mean = mean(dmax)

limiar_p_mean = mean(limiar_p)
dmax_p_mean = mean(dmax_p)

beep

%graficos e analises antigas do HD e HD-PATV, verificar se pode ser util
%dmax_min = min(dmax)
% dmax_mean = mean(dmax)
% dmax_std = std(dmax)
%createfigure1(SNR, HE_eps')
%createfigure2(SNR, HE_eps')