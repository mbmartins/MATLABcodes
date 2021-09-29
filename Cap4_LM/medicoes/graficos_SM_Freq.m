clear all; close all; clc;

load('Dados.mat')

freq_P(:,1) = 1e-6*DadosSMFreq.t1;
freq_P(:,2) = 1e-6*DadosSMFreq.t2;
freq_P(:,3) = 1e-6*DadosSMFreq.t3;
freq_P(:,4) = 1e-6*DadosSMFreq.t4;
freq_P(:,5) = 1e-6*DadosSMFreq.t5;
freq_P(:,6) = 1e-6*DadosSMFreq.t6;
freq_P(:,7) = 1e-6*DadosSMFreq.t7;
freq_P(:,8) = 1e-6*DadosSMFreq.t8;
freq_P(:,9) = 1e-6*DadosSMFreq.t9;

freq_P = freq_P*60; % para ficar em [Hz]

freq_nom = 60;
casos = DadosSMFreq.caso

diff = freq_P;

um = freq_nom*1e-5;
Um = um*ones(1,9);

T = 0.1; %[s]
tau_s = T*(0.1:0.1:0.9); %[s]

f = figure(1)
f.Position = [293 478.6000 754.4000 283.2000];
    errorbar(tau_s,diff(1,:),Um,Um,'bo','LineWidth',1.0,'MarkerSize', 2)
    hold on;
    errorbar(tau_s,diff(2,:),Um,Um,'ro','LineWidth',1.0,'MarkerSize', 2)
    errorbar(tau_s,diff(3,:),Um,Um,'ko','LineWidth',1.0,'MarkerSize', 2)
    difp_mean = (mean([diff(1,:),diff(2,:),diff(3,:)]))*ones(1,length(tau_s));
    difp_std = (std([diff(1,:),diff(2,:),diff(3,:)]))*ones(1,length(tau_s));
    plot(tau_s,difp_mean,'k:','LineWidth',2.0)
    plot(tau_s,difp_mean+difp_std,'m:','LineWidth',2.0)
    plot(tau_s,difp_mean-difp_std,'m:','LineWidth',2.0)
    lg = legend('$\phi_0 = 0^o$','$\phi_0 = 120^o$','$\phi_0 = -120^o$',"$\mu$",'$\mu \pm \sigma$');
    lg.Interpreter = 'latex';
    lg.Location = 'eastoutside';
    lg.Orientation = 'vertical';
    lg.FontSize = 12;
    grid on;
    xlb = xlabel('$\tau$ [s]','Interpreter','latex')
    xlb.FontSize = 14;
    xlim([0 .1]); %ylim([-1e-2 1e-2]);
    ylb = ylabel('$\epsilon_{f} = (\hat{f} - f_{ref})$ [Hz]','Interpreter','latex')
    ylb.FontSize = 14;
%    title('$h_m=0,1$','Interpreter','latex')
% 

f2 = figure(2)
f2.Position = [293 478.6000 754.4000 283.2000];
errorbar(tau_s,diff(4,:),Um,'bo')
hold on;
errorbar(tau_s,diff(5,:),Um,'ro')
errorbar(tau_s,diff(6,:),Um,'ko')
legend(casos(4:6))
grid on;
xlabel('$\tau$ relativo a $T$','Interpreter','latex')
xlim([0 .1])
ylabel('$\epsilon_{f} = (f_{nom} - f_{P})$ [Hz]','Interpreter','latex')