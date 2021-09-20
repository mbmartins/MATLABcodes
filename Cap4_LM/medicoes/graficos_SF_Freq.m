clear all; close all; clc;

load('Dados.mat')

freq_P(:,1) = 1e-6*DadosSFFreq.t1;
freq_P(:,2) = 1e-6*DadosSFFreq.t2;
freq_P(:,3) = 1e-6*DadosSFFreq.t3;
freq_P(:,4) = 1e-6*DadosSFFreq.t4;
freq_P(:,5) = 1e-6*DadosSFFreq.t5;
freq_P(:,6) = 1e-6*DadosSFFreq.t6;
freq_P(:,7) = 1e-6*DadosSFFreq.t7;
freq_P(:,8) = 1e-6*DadosSFFreq.t8;
freq_P(:,9) = 1e-6*DadosSFFreq.t9;

freq_P = freq_P*60; % para ficar em [Hz]

freq_nom = 60;
casos = DadosSFFreq.caso

diff = freq_P;

um = freq_nom*1e-5;
Um = um*ones(1,9);

tau_rel = 0.1:0.1:0.9;

figure(1)
subplot(121)
errorbar(tau_rel,diff(1,:),Um,'bo')
hold on;
errorbar(tau_rel,diff(2,:),Um,'ro')
errorbar(tau_rel,diff(3,:),Um,'ko')
legend(casos(1:3))
grid on;
xlabel('$\tau$ relativo a $T$','Interpreter','latex')
xlim([0 1])
ylabel('$\epsilon_{f} = (f_{nom} - f_{P})$ [Hz]','Interpreter','latex')
subplot(122)
h = histogram(diff)
h.Orientation = 'horizontal'

figure(2)
errorbar(tau_rel,diff(4,:),Um,'bo')
hold on;
errorbar(tau_rel,diff(5,:),Um,'ro')
errorbar(tau_rel,diff(6,:),Um,'ko')
legend(casos(4:6))
grid on;
xlabel('$\tau$ relativo a $T$','Interpreter','latex')
xlim([0 1])
ylabel('$\epsilon_{f} = (f_{nom} - f_{P})$ [Hz]','Interpreter','latex')