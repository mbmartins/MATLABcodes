clear all; close all; clc;

load('Dados.mat')

X_P(:,1) = DadosSMMag.t1;
X_P(:,2) = DadosSMMag.t2;
X_P(:,3) = DadosSMMag.t3;
X_P(:,4) = DadosSMMag.t4;
X_P(:,5) = DadosSMMag.t5;
X_P(:,6) = DadosSMMag.t6;
X_P(:,7) = DadosSMMag.t7;
X_P(:,8) = DadosSMMag.t8;
X_P(:,9) = DadosSMMag.t9;

X_nom = 1;
% casos = DadosSMMag.caso

diff = X_nom - X_P;

um = 152e-6;
Um = um*ones(1,9);

tau_rel = 0.1:0.1:0.9;

figure(1)
errorbar(tau_rel,diff(1,:),Um,'bo')
hold on;
errorbar(tau_rel,diff(2,:),Um,'ro')
errorbar(tau_rel,diff(3,:),Um,'ko')
% legend(casos(1:3))
grid on;
xlabel('$\tau$ relativo a $T$','Interpreter','latex')
xlim([0 1])
ylabel('$\epsilon_{X} = (X_{nom} - X_{P})$ [V]','Interpreter','latex')

figure(2)
errorbar(tau_rel,diff(4,:),Um,'bo')
hold on;
errorbar(tau_rel,diff(5,:),Um,'ro')
errorbar(tau_rel,diff(6,:),Um,'ko')
% legend(casos(4:6))
grid on;
xlabel('$\tau$ relativo a $T$','Interpreter','latex')
xlim([0 1])
ylabel('$\epsilon_{X} = (X_{nom} - X_{P})$ [V]','Interpreter','latex')