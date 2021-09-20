clear all; close all; clc;

load('Dados.mat')

phi_P(:,1) = DadosSMFase.t1;
phi_P(:,2) = DadosSMFase.t2;
phi_P(:,3) = DadosSMFase.t3;
phi_P(:,4) = DadosSMFase.t4;
phi_P(:,5) = DadosSMFase.t5;
phi_P(:,6) = DadosSMFase.t6;
phi_P(:,7) = DadosSMFase.t7;
phi_P(:,8) = DadosSMFase.t8;
phi_P(:,9) = DadosSMFase.t9;

phi_nom = DadosSMFase.Phase;
casos = DadosSMFase.caso

diff = phi_nom - phi_P;

um = 6e-3;
Um = um*ones(1,9);

tau_rel = 0.1:0.1:0.9;

figure(1)
errorbar(tau_rel,diff(1,:),Um,'bo')
hold on;
errorbar(tau_rel,diff(2,:),Um,'ro')
errorbar(tau_rel,diff(3,:),Um,'ko')
legend(casos(1:3))
grid on;
xlabel('$\tau$ relativo a $T$','Interpreter','latex')
xlim([0 1])
ylabel('$\epsilon_{\phi} = (\phi_{nom} - \phi_{P})$ [graus]','Interpreter','latex')

figure(2)
errorbar(tau_rel,diff(4,:),Um,'bo')
hold on;
errorbar(tau_rel,diff(5,:),Um,'ro')
errorbar(tau_rel,diff(6,:),Um,'ko')
legend(casos(4:6))
grid on;
xlabel('$\tau$ relativo a $T$','Interpreter','latex')
xlim([0 1])
ylabel('$\epsilon_{\phi} = (\phi_{nom} - \phi_{P})$ [graus]','Interpreter','latex')