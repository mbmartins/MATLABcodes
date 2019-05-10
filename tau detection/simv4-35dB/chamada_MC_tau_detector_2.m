clear all; close all; clc;

%chamada para MC - detecção de tau 

SNR = 35; %SNR em db
Pin = 90; %fase inicial em graus

%conjuntos de limiares para a simulação
%lim_mag = 0.1:0.1:7; %limiar do detector de magnitude (x mediana)
%lim_fase = 1:0.2:7;%limiar do detector de fase (x mediana)
%obs: nos valores abaixo, praticamente atua somente o detector de mag
lim_mag = 1; %limiar do detector de magnitude (x mediana)
lim_fase = 7;%limiar do detector de fase (x mediana)

%Nruns = 10000; %numero de rodadas de Monte Carlo
Nruns = 1;
KxS = 0.1;

for m = 1:length(lim_mag)
    for n = 1:length(lim_fase)
        %tau_pp_lims = [0.4990 0.5010]; %conjunto de posições relativas de tau dentro da janela
        tau_pp_lims = [0.2 0.2];
        [tau_errors,tau_pp, extremos(m,n), det_mag(m,n),det_fase(m,n),det_nan(m,n)] = MC_tau_error_hibrid_detector_v4_taurand(KxS,SNR,Pin,lim_mag(m),lim_fase(n),tau_pp_lims,Nruns);
        display([lim_mag(m) lim_fase(n)])
        
    end
    tau_vs_errors(m,:) = [tau_pp extremos(m,n)];
end

dt = 1/5000;
plot(tau_pp,tau_errors/dt,'.'); xlabel('\tau [% T]'); ylabel('\epsilon_{\tau} [\Delta t]')
figure
%obs: aplicando a correção sistematica
corr = -1;
histfit(tau_errors/dt+corr,10); xlabel('\epsilon_{\tau} [\Deltat]')
beep
erro_tau_percentual=extremos*100/Nruns;
det_mag_perc = det_mag*100/Nruns;
det_fase_perc = det_fase*100/Nruns;
det_nan_perc = det_nan*100/Nruns;



display('OK')

%surf(erro_tau_percentual)
% 
surf(lim_mag,lim_fase,erro_tau_percentual')
xlabel('lim_{mag}');ylabel('lim_{fase}')