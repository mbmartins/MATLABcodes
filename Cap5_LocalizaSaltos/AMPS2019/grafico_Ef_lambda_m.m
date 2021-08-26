close all; clear all; clc;
%desenha a figura do grafico de E_f em função de lambda

load('otimiza_lambda_f.mat')
figure(1)
%subplot(1,2,1)
plot(lambda_f,eps_HD_PATV); hold on;
plot(lambda_f,4.55*ones(1,length(lambda_f)),'k');
xlabel('$\lambda_f$','Interpreter','latex'); 
%ylabel('$\mathcal{E}_f$ [%]','Interpreter','latex')
ylabel('$\mathcal{E}_f$ [\%]','Interpreter','latex')
legend('\epsilon > 2 \Delta t','Desempenho aceitável')
xlim([min(lambda_f) max(lambda_f)])
grid on;

%desenha a figura do grafico de E_f em função de lambda
% load('otimiza_lambda_m.mat')
% figure(2)
% %subplot(1,2,1)
% plot(lambda_m,eps_HD_PATV); hold on;
% plot(lambda_m,4.55*ones(1,length(lambda_m)),'k');
% xlabel('$\lambda_f$','Interpreter','latex'); 
% %ylabel('$\mathcal{E}_f$ [%]','Interpreter','latex')
% ylabel('$\mathcal{E}_f$ [\%]','Interpreter','latex')
% legend('\epsilon > 2 \Delta t','Desempenho aceitável')
% xlim([min(lambda_m) max(lambda_m)])
% grid on;
