clear all; close all; clc;
load('MC_localiza_salto_freq10000eps8')
h3 = histogram(tau_error_FD, 'FaceAlpha',0.5); hold on;
h4 = histogram(tau_error_FD_PATV, 'FaceAlpha',0.5); hold on;
h3.Normalization = 'probability';h3.BinWidth = 8;
h4.Normalization = 'probability';h4.BinWidth = 8;
h3.EdgeColor = 'none'; h4.EdgeColor = 'none';
ylabel('Frequência relativa','Fontsize',14)
xlabel('\epsilon_{\tau}','Fontsize',14)%,'Interpreter','latex')
grid on
line([-8,-8],[0,0.5])
line([8,8],[0,0.5])
%ylim([-240 240]);
% 
% legend('Interpreter','latex')
legend('EF','EF-PATV');
% legend('Location','northwest');
% legend('EdgeColor','none');
% legend('Color','none');
% legend('Orientation','horizontal');
% legend('Fontsize',12)
