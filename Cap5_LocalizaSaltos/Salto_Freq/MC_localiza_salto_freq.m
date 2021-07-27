%generation of performance comparison
clear all; close all; clc

%OBS: resultados pessimos, não vale a pena insistir muito
% apostar no PATV

SNR = 60;
%SNR = 60:-10:30;
%SNR = 90:-10:30;

%angulo phi_0, fixo ou tabela
Pss = 30;  %pior caso para freq é phi0 = 0;
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

%fator de multiplicação para limiares de detecção
% km = 3;
kf = 10; % para salto de frequencia pq 10x??
% 
% parametro para PATV_HE
lambda = 2.5;

MCruns = 100;

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
             
            [tau_e,dmax(r),limiar(r)] = HF_estimator(Signal,kf); %kf é o limiar de detecção
            [tau_ep,dmax_p(r),limiar_p(r)] = HF_PATV_estimator(Signal,kf,lambda); %kf é o limiar de detecção

            %tau_e é dado em dt
            %tau_e = tau_e/Fs; dt=1/Fs;
            NSamples = length(Signal);
            tau_error_HF(r) = (tau_e - floor(tau*NSamples)); %erro desta rodada de MC
            tau_error_HF_PATV(r) = (tau_ep - floor(tau*NSamples)); %erro desta rodada de MC
            
        %Performance of tau estimation    
        crit_tau_HE(r) = abs(tau_error_HF(r))>2;
        crit_nan_HE(r) = isnan(tau_error_HF(r));
        crit_HE(r) = (crit_tau_HE(r) || crit_nan_HE(r));
        
        %performance para PATV
%          crit_tau_PATV_HE(r) = abs(tau_error_PATV_HE(r,1))>2;
%          crit_nan_PATV_HE(r) = isnan(tau_error_PATV_HE(r,1));
%          crit_PATV_HE(r) = (crit_tau_PATV_HE(r) || crit_nan_PATV_HE(r));
        
        end
        false_pos_HE(ps,s) = sum(crit_tau_HE)*100/MCruns;
        false_neg_HE(ps,s) = sum(crit_nan_HE)*100/MCruns;
%          false_pos_PATV_HE(ps,s) = sum(crit_tau_PATV_HE)*100/MCruns;
%          false_neg_PATV_HE(ps,s) = sum(crit_nan_PATV_HE)*100/MCruns;
        
        eps_HE(ps,s) = sum(crit_HE)*100/MCruns; % in [%];
%         eps_PATV_HE(ps,s) = sum(crit_PATV_HE)*100/MCruns; % in [%]; 
    message = "phi:"+Ps+"°,  SNR:"+SNR(s)+" dB, MCrun:"+r
    end
end

limiar_mean = mean(limiar)
dmax_mean = mean(dmax)
beep

%graficos e analises antigas do HD e HD-PATV, verificar se pode ser util
%dmax_min = min(dmax)
% dmax_mean = mean(dmax)
% dmax_std = std(dmax)
%createfigure1(SNR, HE_eps')
%createfigure2(SNR, HE_eps')