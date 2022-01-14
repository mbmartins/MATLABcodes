close all; clear all; clc;

load('salto_freq_lambdaMC_MedSF300_v2.mat')

FEmax5 = max(abs(FE5'));
FE_std5 = std(FE5');
FEmean5 = abs(mean(FE5'));

FE1max5 = max(abs(FE15'));
FE1_std5 = std(FE15');
FE1mean5 = abs(mean(FE15'));

KFEmax5 = max(abs(kfE5'));
KFE_std5 = std(kfE5');
KFEmean5 = abs(mean(kfE5'));

close all;
%----- Figuras para FE x lambda

figs(2) = figure(2);hold off

%subplot(131)
semilogy(lambda_n,abs(FEmean5),'b.-');hold on; 
semilogy(lambda_n,abs(FEmean5) + FE_std5,'b--');
semilogy(lambda_n,abs(FE1mean5),'r.-'); hold on;
semilogy(lambda_n,abs(FE1mean5)+FE1_std5,'r--');
semilogy(lambda_n,abs(KFEmean5),'k.-'); hold on;
semilogy(lambda_n,abs(KFEmean5)+KFE_std5,'k--');

xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
xlabel('$\lambda$','Fontsize',16)
ylabel('$|FE|,|FE_1|,|h_fE|$ [Hz]','Fontsize',15)
legend('|\mu(FE)|','|\mu(FE)| + \sigma(FE)',...
       '|\mu(FE_1)|','|\mu(FE_1)| + \sigma(FE_1)',...
       '|\mu(h_fE)|','|\mu(h_fE)| + \sigma(h_fE)')
legend('Location','eastoutside')
legend('Fontsize',12)
%xlabel('Fontsize',14)
%ylabel('Fontsize',14)
legend('Orientation','vertical')
grid on
ylim([10e-3 1]);
xlim([0.15 3])
