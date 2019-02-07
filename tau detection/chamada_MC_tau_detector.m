clear all; close all; clc;

%chamada para MC - detecção de tau 

SNR = 60; %SNR em db
Pin = 90; %fase inicial em graus
lim_mag = 1:7; %limiar do detector de magnitude (x mediana)
lim_fase = 5:11;%limiar do detector de fase (x mediana)
tau_pp_set = 0.5; %conjunto de posições relativas de tau dentro da janela
Nruns = 10000; %numero de rodadas de Monte Carlo

for m = 1:length(lim_mag)
    for n = 1:length(lim_fase)
        [tau_errors, extremos(m,n), det_mag(m,n),det_fase(m,n),det_nan(m,n)] = MC_tau_error_hibrid_detector(SNR,Pin,lim_mag(m),lim_fase(n),tau_pp_set,Nruns);
    end
end

beep
erro_tau_percentual=extremos*100/Nruns;
det_mag_perc = det_mag*100/Nruns;
det_fase_perc = det_fase*100/Nruns;
det_nan_perc = det_nan*100/Nruns;

display('OK')
% 
% surf(lim_mag,lim_fase,erro_tau_percentual')
% xlabel('lim_{mag}');ylabel('lim_{fase}')