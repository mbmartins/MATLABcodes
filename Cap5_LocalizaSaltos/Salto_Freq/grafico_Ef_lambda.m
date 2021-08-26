close all; clear all; clc;
%desenha a figura do grafico de E_f em função de lambda

load('otimiza_lambda_r.mat')
figure(1)
subplot(1,2,1)
semilogy(lambda,eps_FD_PATV); hold on;
semilogy(lambda,eps_FD_PATV_7);
semilogy(lambda,eps_FD_PATV_8);
semilogy(lambda,4.55*ones(1,length(lambda)),'k');
xlabel('$\lambda_r$','Interpreter','latex'); 
%ylabel('$\mathcal{E}_f$ [%]','Interpreter','latex')
ylabel('$\mathcal{E}_f$ [\%]','Interpreter','latex')
legend('\epsilon > 2 \Delta t','\epsilon > 7 \Delta t','\epsilon > 8 \Delta t','Desempenho aceitável')
xlim([0.02 max(lambda)])
grid on;
subplot(1,2,2)
semilogy(lambda,std_tau_e_FD_PATV);
xlabel('\lambda_r'); ylabel('\sigma_{\tau}')
xlim([0.02 max(lambda)])
grid on