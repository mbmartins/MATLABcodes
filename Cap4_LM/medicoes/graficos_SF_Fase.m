clear all; close all; clc;

load('Dados.mat')

phi_P(:,1) = DadosSFFase.t1(5);
phi_P(:,2) = DadosSFFase.t2(5);
phi_P(:,3) = DadosSFFase.t3(5);
phi_P(:,4) = DadosSFFase.t4(5);
phi_P(:,5) = DadosSFFase.t5(5);
phi_P(:,6) = DadosSFFase.t6(5);
phi_P(:,7) = DadosSFFase.t7(5);
phi_P(:,8) = DadosSFFase.t8(5);
phi_P(:,9) = DadosSFFase.t9(5);

tau_rel = 0.1:0.1:0.9;

ha = 10;
phi_1 = 120;
phi_2 = 120 + ha;

phi_ref = phi_1*tau_rel + phi_2*(1-tau_rel);

diff = phi_P - phi_ref;

um = 15e-3; % 15 m graus
Um = um*ones(1,9);
T = 0.1;
tau_rel = 0.1:0.1:0.9;
tau_s = tau_rel*T;

f = figure(1)
f.Position = [293 478.6000 754.4000 283.2000];
    errorbar(tau_s,diff(1,:),Um,Um,'bo','LineWidth',1.0,'MarkerSize', 2)
    hold on;
%     errorbar(tau_s,diff(2,:),Um,Um,'ro','LineWidth',1.0,'MarkerSize', 2)
%     errorbar(tau_s,diff(3,:),Um,Um,'ko','LineWidth',1.0,'MarkerSize', 2)
    difp_mean = (mean([diff(1,:)]))*ones(1,length(tau_s));
    difp_std = (std([diff(1,:)]))*ones(1,length(tau_s));
    plot(tau_s,difp_mean,'k:','LineWidth',2.0,'MarkerSize', 2)
    plot(tau_s,difp_mean+difp_std,'m:','LineWidth',2.0,'MarkerSize', 2)
    plot(tau_s,difp_mean-difp_std,'m:','LineWidth',2.0,'MarkerSize', 2)
    lg = legend('$\phi_0 = 120^o$',"$\mu$",'$\mu \pm \sigma$');
    lg.Interpreter = 'latex';
    lg.Location = 'eastoutside';
    lg.Orientation = 'vertical';
    lg.FontSize = 12;
    grid on;
    xlb = xlabel('$\tau$ [s]','Interpreter','latex')
    xlb.FontSize = 14;
    xlim([0 .1]); %ylim([-1e-2 1e-2]);
    ylb = ylabel('$\epsilon_{\phi_0} = (\hat{\phi_0} - \phi_0)$ [graus]','Interpreter','latex')
    ylb.FontSize = 14;
% 
% figure(2)
% errorbar(tau_rel,diff(4,:),Um,'bo')
% hold on;
% errorbar(tau_rel,diff(5,:),Um,'ro')
% errorbar(tau_rel,diff(6,:),Um,'ko')
% legend(casos(4:6))
% grid on;
% xlabel('$\tau$ relativo a $T$','Interpreter','latex')
% xlim([0 1])
% ylabel('$\epsilon_{\phi} = (\phi_{nom} - \phi_{P})$ [graus]','Interpreter','latex')