close all; clear all; clc;
%desenha a figura do grafico de E_f em função de lambda

load('otimiza_lambda_r2.mat')
%ajuste de escalas
lambda = 0.02:0.02:0.4;
indices = 1:2:length(eps_FD_PATV);
eps_FD_PATV_down = eps_FD_PATV(indices);
eps_FD_PATV_7_down = eps_FD_PATV_7(indices);
eps_FD_PATV_8_down = eps_FD_PATV_8(indices);

figure(1)
%subplot(1,2,1)
semilogy(lambda,eps_FD_PATV_down); hold on;
semilogy(lambda,eps_FD_PATV_7_down);
semilogy(lambda,eps_FD_PATV_8_down);
semilogy(lambda,4.55*ones(1,length(lambda)),'k');
xlabel('$\lambda_r$','Interpreter','latex'); 
%ylabel('$\mathcal{E}_f$ [%]','Interpreter','latex')
ylabel('$\mathcal{E}_f$ [\%]','Interpreter','latex')
legend('\epsilon > 2 \Delta t','\epsilon > 7 \Delta t','\epsilon > 8 \Delta t','Desempenho aceitável')
xlim([0.02 max(lambda)])
grid on;
% subplot(1,2,2)
% semilogy(lambda,std_tau_e_FD_PATV);
% xlabel('\lambda_r'); ylabel('\sigma_{\tau}')
% xlim([0.02 max(lambda)])
% grid on