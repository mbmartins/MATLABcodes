clear all; close all; clc;

load('Dados.mat')
clear X_P;

X_P(:,1) = DadosSFMag.t1;
X_P(:,2) = DadosSFMag.t2;
X_P(:,3) = DadosSFMag.t3;
X_P(:,4) = DadosSFMag.t4;
X_P(:,5) = DadosSFMag.t5;
X_P(:,6) = DadosSFMag.t6;
X_P(:,7) = DadosSFMag.t7;
X_P(:,8) = DadosSFMag.t8;
X_P(:,9) = DadosSFMag.t9;

tau_rel = 0.1:0.1:0.9;

hm = 0;
X1 = 1;
X2 = 1 + hm;
X_nom = (X1*tau_rel + (1-tau_rel)*X2);
% casos = DadosSFMag.caso

diff = X_P - X_nom;

um = 152e-6;
Um = um*ones(1,9);

T = 0.1;
tau_s = tau_rel*T; %segundos

diff_mean = mean([diff(:,1);diff(:,2);diff(:,3)])*ones(1,length(tau_s));
diff_std = std([diff(:,1);diff(:,2);diff(:,3)])*ones(1,length(tau_s));

f = figure(1)
f.Position = [293 478.6000 754.4000 283.2000];
errorbar(tau_s,diff(1,:),Um,'bo','LineWidth',1.0,'MarkerSize', 2)
hold on;
errorbar(tau_s,diff(2,:),Um,'ro','LineWidth',1.0,'MarkerSize', 2)
errorbar(tau_s,diff(3,:),Um,'ko','LineWidth',1.0,'MarkerSize', 2)
plot(tau_s,diff_mean,'k:','LineWidth',2.0)
plot(tau_s,diff_mean + diff_std,'m:','LineWidth',2.0)
plot(tau_s,diff_mean - diff_std,'m:','LineWidth',2.0)
% errorbar(tau_rel,diff(2,:),Um,'ro')
% errorbar(tau_rel,diff(3,:),Um,'ko')
% legend(casos(1:3))
grid on;
xlb = xlabel('$\tau$ [s]','Interpreter','latex')
xlim([0 .1])
    xlb.FontSize = 14;
ylb = ylabel('$\epsilon_{|X|} = (|\hat{X}| - |X_{ref|})$ [V]','Interpreter','latex')
ylb.FontSize = 14;
ylim([-5e-3 -3.2e-3])

lg = legend('$\phi_0 = 0^o$','$\phi_0 = 120^o$','$\phi_0 = -120^o$',"$\mu$",'$\mu \pm \sigma$');
    lg.Interpreter = 'latex';
    lg.Location = 'eastoutside';
    lg.Orientation = 'vertical';
    lg.FontSize = 12;