clear all; close all; clc;

load('Dados.mat')

X_P(1) = DadosSMMag.t1(5);
X_P(2) = DadosSMMag.t2(5);
X_P(3) = DadosSMMag.t3(5);
X_P(4) = DadosSMMag.t4(5);
X_P(5) = DadosSMMag.t5(5);
X_P(6) = DadosSMMag.t6(5);
X_P(7) = DadosSMMag.t7(5);
X_P(8) = DadosSMMag.t8(5);
X_P(9) = DadosSMMag.t9(5);

tau_rel = 0.1:0.1:0.9;

hm = -0.1;
X1 = 1;
X2 = 1 + hm;
X_nom = (X1*tau_rel + (1-tau_rel)*X2);
% casos = DadosSMMag.caso

diff = X_P - X_nom;

um = 152e-6;
Um = um*ones(1,9);

T = 0.1;
tau_s = tau_rel*T; %segundos

diff_mean = mean(diff)*ones(1,length(tau_s));
diff_std = std(diff)*ones(1,length(tau_s));

f = figure(1)
f.Position = [293 478.6000 754.4000 283.2000];
errorbar(tau_s,diff(1,:),Um,'bo','LineWidth',1.0,'MarkerSize', 2)
hold on;
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

lg = legend('$\phi_0 = 0^o$',"$\mu$",'$\mu \pm \sigma$');
    lg.Interpreter = 'latex';
    lg.Location = 'eastoutside';
    lg.Orientation = 'vertical';
    lg.FontSize = 12;