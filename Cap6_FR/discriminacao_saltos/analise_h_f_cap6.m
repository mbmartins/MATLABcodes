clear all; close all; clc;
load('boxplot_ROCOF_geral_kfest.mat')

duracao = (tau_vec(2,1) - tau_vec(1,1))*F0*T; % em ciclos de F0

tau_inicial = tau_vec(1,:); % de 0.1 a 0.9?? não faz sentido...

%melhor seria limitar a duracao a 1 ciclo 
% o tau inicial começa em t1 = 0.1 ok
% mas vai até t1 = 0.9 - 1 ciclo

%define o discriminador para h_f > 0.5

% dá para aproveitar os dados ignorando os pontos cujos t1 > 0.9 - 1.8
% ciclo

i = 5
limiar = 0.5;
lim_vec = limiar*ones(1,length(tau_vec));

f2 = figure(2);
% subplot(121)
f2.Position = [ 143.4000  419.4000  624.0000  342.4000]
plot(tau_vec(1,:)*T,lim_vec,'r'); hold on;
plot(tau_vec(1,:)*T,abs(kfest_mag(i,:)),'b.')
plot(tau_vec(1,:)*T,abs(kfest_fase(i,:)),'m.')
hold on; 
plot(tau_vec(1,:)*T,abs(kfest_freq(i,:)),'k.')
%p.Interpreter = 'latex'
%p.Interpreter = 'latex'
xlb = xlabel('$\tau$ [s]'); xlb.Interpreter = 'latex'; xlb.FontSize = 14;
ylb = ylabel('$|\hat{h}_f|$ [Hz]'); ylb.Interpreter = 'latex'; ylb.FontSize = 14;
lgd = legend('Limiar detecção','Salto Magnitude','Salto Fase','Salto Frequência');
lgd.FontSize = 12;
lgd.Location = 'north'
grid on;

% subplot(122)
% f2.Position = [143.4000  324.2000  973.6000  437.6000]
% plot(tau_vec(1,:)*T,lim_vec,'r'); hold on;
% plot(tau_vec(1,:)*T,abs(kfest_fase(i,:)),'b.')
% hold on; 
% plot(tau_vec(1,:)*T,abs(kfest_freq(i,:)),'k.')
% xlb = xlabel('$\tau$ [s]'); xlb.Interpreter = 'latex'; xlb.FontSize = 14;
% ylb = ylabel('$|\hat{h}_f|$ [Hz]'); ylb.Interpreter = 'latex'; ylb.FontSize = 14;
% lgd = legend('Limiar detecção','Salto Fase','Salto Frequência');
% lgd.FontSize = 14;
% lgd.Location = 'north'
% grid on;

%title("h_f estimado por MedSF-PATV")
print(f2,'hf_discrim','-dpng','-r600')

limiar= 0.5;

% quantos pontos de salto de frequencia h_f < limiar?

Falhas_freq = 100*sum(abs(kfest_freq(i,:))<limiar)/MCiter  % em porcento

% quantos pontos de afundamento teve h_f > limiar?
Falhas_sag = 100*sum(abs(kfest_sag(i,:))>limiar)/MCiter  % em porcento

print(f2,'hf_disc','-dpng','-r600')

Falhas_freq_restrito = 100*sum(abs(kfest_freq(i,:))<limiar)/MCiter  % em porcento

% quantos pontos de afundamento teve h_f > limiar?
Falhas_sag_restrito = 100*sum(abs(kfest_sag(i,:))>limiar)/MCiter  % em porcento