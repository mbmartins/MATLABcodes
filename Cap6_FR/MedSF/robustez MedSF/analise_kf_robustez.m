clear all; close all; clc;

%dados provenientes de estimadores_ROCOF_geral.m
load('boxplot_ROCOF_geral_kfest.mat')

figure(1)
subplot(2,1,1)
%MedSF = EF1
h3 = histogram(kfest_fase(1,:), 'FaceAlpha',0.5); hold on;
h4 = histogram(kfest_freq(5,:), 'FaceAlpha',0.5); hold on;
h3.Normalization = 'probability';h2.BinWidth = 0.08;
h4.Normalization = 'probability';h4.BinWidth = 0.08;
h3.EdgeColor = 'none'; h4.EdgeColor = 'none';
ylabel('Frequ�ncia relativa','Fontsize',14)
ylim([0 0.5]);
%xlabel('Interpreter','latex')
% t = xlabel('$\hat{h}_f$ [Hz]','Fontsize',14);
% t.Interpreter = 'latex';
%xlim([-0.5 1.5])
legend('Interpreter','latex')
legend('Salto Magnitude ($h_m = -0,1$)','Salto Frequ\^encia ($h_f = 1$ Hz)');
legend('Location','northwest');
legend('EdgeColor','none');
legend('Color','none');
legend('Orientation','horizontal');
legend('Fontsize',12)
title('Salto de Magnitude - MedSF')

subplot(2,1,2)
%MedSF-PATV = EF5
h3 = histogram(kfest_fase(5,:), 'FaceAlpha',0.5); hold on;
h4 = histogram(kfest_freq(5,:), 'FaceAlpha',0.5); hold on;
h3.Normalization = 'probability';h2.BinWidth = 0.08;
h4.Normalization = 'probability';h4.BinWidth = 0.08;
h3.EdgeColor = 'none'; h4.EdgeColor = 'none';
ylabel('Frequ�ncia relativa','Fontsize',14);
ylim([0 0.5]);
t = xlabel('$\hat{h}_f$ [Hz]','Fontsize',14);
t.Interpreter = 'latex';
%xlim([-0.5 1.5])
legend('Interpreter','latex')
legend('Salto Magnitude ($h_m = -0,1$)','Salto Frequ\^encia ($h_f = 1$ Hz)');
legend('Location','northwest');
legend('EdgeColor','none');
legend('Color','none');
legend('Orientation','horizontal');
legend('Fontsize',12)
title('Salto de Magnitude - MedSF-PATV')

figure(2)
subplot(2,1,1)
h2 = histogram(kfest_mag(1,:), 'FaceAlpha',0.5); hold on;
h4 = histogram(kfest_freq(5,:), 'FaceAlpha',0.5); hold on;
h2.Normalization = 'probability';h2.BinWidth = 0.08;
h4.Normalization = 'probability';h4.BinWidth = 0.08;
h2.EdgeColor = 'none'; h4.EdgeColor = 'none';
%xlim([-0.5 1.5])
ylabel('Frequ�ncia relativa','Fontsize',14);
ylim([0 0.5]);
% t = xlabel('$\hat{h}_f$ [Hz]','Fontsize',14);
% t.Interpreter = 'latex';
legend('Interpreter','latex')
legend('Salto Fase ($h_a = 10^o$)','Salto Frequ\^encia ($h_f = 1$ Hz)');
legend('Location','northeast');
legend('Orientation','horizontal');
legend('EdgeColor','none');
legend('Color','none');
legend('Fontsize',12)
title('Salto de Fase - MedSF')

subplot(2,1,2)
h2 = histogram(kfest_mag(5,:), 'FaceAlpha',0.5); hold on;
h4 = histogram(kfest_freq(5,:), 'FaceAlpha',0.5); hold on;
h2.Normalization = 'probability';h2.BinWidth = 0.08;
h4.Normalization = 'probability';h4.BinWidth = 0.08;
h2.EdgeColor = 'none'; h4.EdgeColor = 'none';
%xlim([-0.5 1.5])
ylabel('Frequ�ncia relativa','Fontsize',14);
ylim([0 0.5]);
t = xlabel('$\hat{h}_f$ [Hz]','Fontsize',14);
t.Interpreter = 'latex';
legend('Interpreter','latex')
legend('Salto Fase ($h_a = 10^o$)','Salto Frequ\^encia ($h_f = 1$ Hz)');
legend('Location','northeast');
legend('Orientation','horizontal');
legend('EdgeColor','none');
legend('Color','none');
legend('Fontsize',12)
title('Salto de Fase - MedSF-PATV')

% figure(3)
% %tentar fazer numa figura s�?
% h2 = histogram(kfest_mag(5,:), 'FaceAlpha',0.5); hold on;
% h3 = histogram(kfest_fase(5,:), 'FaceAlpha',0.5); hold on;
% h4 = histogram(kfest_freq(5,:), 'FaceAlpha',0.5); hold on;
% h2.Normalization = 'probability';h2.BinWidth = 0.08;
% h3.Normalization = 'probability';h3.BinWidth = 0.08;
% h4.Normalization = 'probability';h4.BinWidth = 0.08;
% h2.EdgeColor = 'none'; 
% h3.EdgeColor = 'none'; 
% h4.EdgeColor = 'none';
% %xlim([-0.5 1.5])
% ylabel('Frequ�ncia relativa','Fontsize',14);
% ylim([0 0.5]);
% t = xlabel('$\hat{h}_f$ [Hz]','Fontsize',14);
% t.Interpreter = 'latex';
% legend('Interpreter','latex')
% legend('Salto Mag','Salto Fase ($h_a = 10^o$)','Salto Frequ\^encia ($h_f = 1$ Hz)');
% legend('Location','northeast');
% legend('Orientation','horizontal');
% legend('EdgeColor','none');
% legend('Color','none');
% legend('Fontsize',12)
% title('MedSF-PATV')



%analise dos desvios padr�o
hf_std_MedSF_fase = std(kfest_fase(1,:))
hf_std_MedSF_PATV_fase = std(kfest_fase(5,:))
hf_std_MedSF_mag = std(kfest_mag(1,:))
hf_std_MedSF_PATV_mag = std(kfest_mag(5,:))
hf_std_MedSF_fase = std(kfest_fase(1,:))
hf_std_MedSF_PATV_fase = std(kfest_fase(5,:))
hf_std_MedSF_freq = std(kfest_freq(1,:))
hf_std_MedSF_PATV_freq = std(kfest_freq(5,:))
