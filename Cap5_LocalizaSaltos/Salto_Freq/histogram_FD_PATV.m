clear all; close all; clc;
load('MC_localiza_salto_freq10000eps8')
h3 = histogram(tau_error_FD, 'FaceAlpha',0.5); hold on;
h4 = histogram(tau_error_FD_PATV, 'FaceAlpha',0.5); hold on;
h3.Normalization = 'probability';h2.BinWidth = 2;
h4.Normalization = 'probability';h4.BinWidth = 2;
h3.EdgeColor = 'none'; h4.EdgeColor = 'none';
ylabel('Frequência relativa','Fontsize',14)
xlabel('\epsilon_{\tau}')
%ylim([-240 240]);
% 
% legend('Interpreter','latex')
legend('FD','FD-PATV');
% legend('Location','northwest');
% legend('EdgeColor','none');
% legend('Color','none');
% legend('Orientation','horizontal');
% legend('Fontsize',12)
