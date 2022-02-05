close all; clear all; clc;

load('consenso_MedSFPATV_sem_comp_low.mat')

close all;
%----- Figuras para FE x lambda

plot(lambda_n,abs(FEmean),'b.-');hold on; 
%plot(lambda_n,abs(FEmean) + FE_std,'b-.');
plot(lambda_n,abs(FE1mean),'r.-'); hold on;
%plot(lambda_n,abs(FE1mean)+FE1_std,'r-.');
plot(lambda_n,abs(FE2mean),'g.-'); hold on;
%plot(lambda_n,abs(FE2mean)+FE2_std,'g-.');
plot(lambda_n,abs(HFEmean),'k.-'); hold on;
%plot(lambda_n,abs(HFEmean)+HFE_std,'k-.');

xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
xlabel('$\lambda$','Fontsize',16)
ylabel('[Hz]','Fontsize',15)
% legend('|\mu(FE)|',...
%        % '|\mu(FE)| + \sigma(FE)',...
%        '|\mu(FE_1)|', ...
%        %'|\mu(FE_1)| + \sigma(FE_1)',...
%        '|\mu(FE_2)|', ...
%        %'|\mu(FE_2)| + \sigma(FE_2)',...
%        '|\mu(h_fE)|', ...
%        %'|\mu(h_fE)| + \sigma(h_fE)')
legend('|\mu(FE)|','|\mu(FE_1)|','|\mu(FE_2)|','|\mu(h_fE)|')
legend('Location','eastoutside')
legend('Fontsize',12)
%xlabel('Fontsize',14)
%ylabel('Fontsize',14)
legend('Orientation','vertical')
grid on
%ylim([10e-5 1]);
%xlim([0.1 7])
