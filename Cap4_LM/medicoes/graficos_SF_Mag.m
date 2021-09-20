clear all; close all; clc;

load('Dados.mat')

X_P(:,1) = DadosSFMag.t1;
X_P(:,2) = DadosSFMag.t2;
X_P(:,3) = DadosSFMag.t3;
X_P(:,4) = DadosSFMag.t4;
X_P(:,5) = DadosSFMag.t5;
X_P(:,6) = DadosSFMag.t6;
X_P(:,7) = DadosSFMag.t7;
X_P(:,8) = DadosSFMag.t8;
X_P(:,9) = DadosSFMag.t9;

X_nom = 1;

diff = X_nom - X_P(5,:);

um = 50e-6;
Um = um*ones(1,9);

tau_rel = 0.1:0.1:0.9;

figure(1)
subplot(121)
errorbar(tau_rel,diff(1,:),Um,'bo')
hold on;
% errorbar(tau_rel,diff(2,:),Um,'ro')
% errorbar(tau_rel,diff(3,:),Um,'ko')
% legend(casos(1:3))
grid on;
xlabel('$\tau$ relativo a $T$','Interpreter','latex')
xlim([0 1])
ylabel('$\epsilon_{X} = (X_{nom} - X_{P})$ [V]','Interpreter','latex')
subplot(122)
h = histogram(diff)
h.Orientation = 'horizontal'