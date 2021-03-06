clear all; close all; clc;

%chamada para MC - detec��o de tau 

SNR = 50; %SNR em db
Pin = 90; %fase inicial em graus
%lim_mag = 0.1:0.1:7; %limiar do detector de magnitude (x mediana)
%lim_fase = 1:0.2:7;%limiar do detector de fase (x mediana)

%obs: nestes valores, praticamente atua somente o detector de mag
lim_mag = 1; %limiar do detector de magnitude (x mediana)
lim_fase = 7;%limiar do detector de fase (x mediana)

Nruns = 10000; %numero de rodadas de Monte Carlo

for m = 1:length(lim_mag)
    for n = 1:length(lim_fase)
        tau_pp_set = 0.5; %conjunto de posi��es relativas de tau dentro da janela
        [tau_errors, extremos(m,n), det_mag(m,n),det_fase(m,n),det_nan(m,n)] = MC_tau_error_hibrid_detector_v4(SNR,Pin,lim_mag(m),lim_fase(n),tau_pp_set,Nruns);
        display([lim_mag(m) lim_fase(n)])
    end
end

beep
erro_tau_percentual=extremos*100/Nruns;
det_mag_perc = det_mag*100/Nruns;
det_fase_perc = det_fase*100/Nruns;
det_nan_perc = det_nan*100/Nruns;

display('OK')

surf(erro_tau_percentual)
% 
% surf(lim_mag,lim_fase,erro_tau_percentual')
% xlabel('lim_{mag}');ylabel('lim_{fase}')