close all; clear all; clc;
%load('salto_fase_lambdaMC.mat')
%load('salto_mag_lambdaMC.mat')

load('salto_freq_lambdaMC_MedSF300.mat')
%load('salto_freq_lambdaMC.mat')


%FEmax4 = max(abs(FE4'));
FEmax5 = max(abs(FE5'));
%FEmax6 = max(abs(FE6'));
%FE_std4 = std(FE4');
FE_std5 = std(FE5');
%FE_std6 = std(FE6');
%FEmean4 = abs(mean(FE4'));
FEmean5 = abs(mean(FE5'));
%FEmean6 = abs(mean(FE6'));

%FE1max4 = max(abs(FE14'));
FE1max5 = max(abs(FE15'));
%FE1max6 = max(abs(FE16'));
%FE1_std4 = std(FE14');
FE1_std5 = std(FE15');
%FE1_std6 = std(FE16');
%FE1mean4 = abs(mean(FE14'));
FE1mean5 = abs(mean(FE15'));
%FE1mean6 = abs(mean(FE16'));

%KFEmax4 = max(abs(kfE4'));
KFEmax5 = max(abs(kfE5'));
%KFEmax6 = max(abs(kfE6'));
%KFE_std4 = std(kfE4');
KFE_std5 = std(kfE5');
%KFE_std6 = std(kfE6');
%KFEmean4 = abs(mean(kfE4'));
KFEmean5 = abs(mean(kfE5'));
%KFEmean6 = abs(mean(kfE6'));
close all;
%----- Figuras para FE x lambda

% ----- Figuras para EF4
% figs(1) = figure(1);hold off
% subplot(211)
% semilogy(lambda_n,abs(FEmean4),'b.-'); hold on;
% title('EF4 - mediana')
% grid on;
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$|FE|$ [Hz]')
% legend('|\mu(FE)|')
% %title('EF4')
% subplot(212)
% semilogy(lambda_n,FE_std4,'k--');
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$|FE|$ [Hz]')
% grid on;
% legend('\sigma(FE)')

% semilogy(lambda_n,abs(KFEmean4),'b.-'); hold on;
% semilogy(lambda_n,abs(KFEmean4)+KFE_std4,'k--');
% semilogy(lambda_n,(KFEmax4),'k:');
% grid on;
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('EF4 $|K_fE|$ [Hz]')
% legend('|\mu(K_fE)|','|\mu(K_fE)| + \sigma(K_fE)', '|K_fE_{max}|')
%title('EF4')

%--------- Figuras para EF5

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
%title('a)')
% 
% %subplot(132)
% %semilogy(lambda_n,FE1_std5,'k--');
% semilogy(lambda_n,abs(FE1mean5),'b.-'); hold on;
% semilogy(lambda_n,abs(FE1mean5)+FE1_std5,'k--');
% %semilogy(lambda_n,(FE1max5),'k:');
% grid on;
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$|FE_1|$ [Hz]')
% legend('|\mu(FE_1)|','|\mu(FE_1)| + \sigma(FE_1)')%, '|FE_1_{max}|')
% legend('Location','southoutside')
% legend('Orientation','horizontal')
% title('b)')
% 
% subplot(133)
% semilogy(lambda_n,abs(KFEmean5),'b.-'); hold on;
% semilogy(lambda_n,abs(KFEmean5)+KFE_std5,'k--');
% %semilogy(lambda_n,(KFEmax5),'k:');
% grid on;
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$|h_fE|$ [Hz]')
% legend('|\mu(h_fE)|','|\mu(h_fE)| + \sigma(h_fE)')%, '|h_fE_{max}|')
% legend('Location','southoutside')
% legend('Orientation','horizontal')
% title('c)')

%--------- Figuras para EF6
% figs(3) = figure(3); hold off
% %subplot(211)
% semilogy(lambda_n,abs(FE1mean6),'b.-');hold on;
% semilogy(lambda_n,abs(FE1mean6) + FE1_std6,'k--');
% %semilogy(lambda_n,abs(FE1mean6) - FE1_std6,'k--');
% legend('|\mu(FE_1)|','|\mu(FE_1)| + \sigma(FE_1)')
% %legend('|\mu(FE_1)|')
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$|FE_1|$ [Hz]')
% grid on
% title('a)')
% subplot(312)
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$|FE_1|$ [Hz]')
% legend('\sigma(FE_1)')
% grid on
% title('b)')
% % subplot(212)
% semilogy(lambda_n,(FE1max6),'k:');
% grid on;
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$|FE_1|$ [Hz]')
% legend('|FE_1_{max}|')
% title('c)')
% ---- Figuras para FE1 x lambda
% figure
% semilogy(lambda_n,abs(FE1mean6),'o-')
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$FE1_{\mu}$ [Hz]')
% 
% 
% figure
% semilogy(lambda_n,FE2max,'o-')
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$FE2_{max}$ [Hz]')
% legend('\tau = 0.1T','\tau = 0.2T', '\tau = 0.3T','\tau = 0.4T','\tau = 0.5T')