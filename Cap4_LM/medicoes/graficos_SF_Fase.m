clear all; close all; clc;

load('Dados.mat')

k = DadosSFFase.Timestamp;
phi_P(:,1) = DadosSFFase.t1;
phi_P(:,2) = DadosSFFase.t2;
phi_P(:,3) = DadosSFFase.t3;
phi_P(:,4) = DadosSFFase.t4;
phi_P(:,5) = DadosSFFase.t5;
phi_P(:,6) = DadosSFFase.t6;
phi_P(:,7) = DadosSFFase.t7;
phi_P(:,8) = DadosSFFase.t8;
phi_P(:,9) = DadosSFFase.t9;

phi_nom = 120 + 10*(0.9:-0.1:0.1);
phiLM = phi_P(5,:);

diff = phi_nom - phiLM;

um = 15e-3;
Um = um*ones(1,length(phiLM));

tau_rel = 0.1:0.1:0.9;

errorbar(tau_rel,diff,Um,'o')
grid on;
xlabel('$\tau$ relativo a $T$','Interpreter','latex')
xlim([0 1])
ylabel('$\epsilon_{\phi_0} = \phi_{nom} - \phi_{P}$ [graus]','Interpreter','latex')