% influencia de tau
load('estimadoresFR_salto_freq_tau_n.mat');
close all;

c = ['k','b','c','r','m','g'];
tits = ["a)EF1","b)EF2","c)EF3","d)EF4","e)EF5","f)EF6"];

% ------------- FE com ruído de 60 dB -------------
fig1 = figure('Units','normalized','Position',[0 0 0.5 1])
for k = 1:6
subplot(3,2,k)
xlim([0,180]); ylim([-0.06 0.05])
EF(k) = plot(tau_n,FE(:,k),c(k)); hold on;
shade(tau_n,FE(:,k) - FE_std(:,k),[c(k),'-.'],tau_n,FE(:,k) + FE_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
%plot(tau_n,FE(:,k) - FE_std(:,k),[c(k),'--']); 
grid on;
xlim([48,432]); ylim([-0.5 0.2])
% legend(EF, 'EF1','EF2','EF3','EF4','EF5','EF6')
% legend('Orientation','horizontal')
% legend('Location','northoutside')
%legend('boxoff')
%title('FE com ruido SNR = 60 dB')
xlabel('\tau_n [amostras]')
ylabel('FE [Hz]')
title(tits(k));
end
% --------------------------------------------------
savefig(fig1,'salto_freq_FE_tau_n.fig')
saveas(fig1,'salto_freq_FE_tau_n.png')

fig2 = figure('Units','normalized','Position',[0 0 0.5 1])
% ------------- FE1 com ruído de 60 dB -------------
for k = 1:6
subplot(3,2,k)
EF(k) = plot(tau_n,FE1(:,k),c(k)); hold on;
shade(tau_n,FE1(:,k) - f1_std(:,k),[c(k),'-.'],tau_n,FE1(:,k) + f1_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
%plot(tau_n,FE(:,k) - FE_std(:,k),[c(k),'--']); 
grid on;
xlim([48,432]); ylim([-2.5 1.])
% legend(EF, 'EF1','EF2','EF3','EF4','EF5','EF6')
% legend('Orientation','horizontal')
% legend('Location','northoutside')
%legend('boxoff')
%title('FE com ruido SNR = 60 dB')
xlabel('\tau_n [amostras]')
ylabel('FE1 [Hz]')
title(tits(k));
end
% --------------------------------------------------
%subplot(3,2,4); ylim([0.4 0.5])
savefig(fig2,'salto_freq_F1_tau_n.fig')
saveas(fig2,'salto_freq_F1_tau_n.png')

fig3 = figure('Units','normalized','Position',[0 0 0.5 1]);
% ------------- KFE1 com ruído de 60 dB -------------
for k = 1:6
subplot(3,2,k)
EF(k) = plot(tau_n,KFE1(:,k),c(k)); hold on;
shade(tau_n,KFE1(:,k) - kf_std(:,k),[c(k),'-.'],tau_n,KFE1(:,k) + kf_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
%plot(tau_n,FE(:,k) - FE_std(:,k),[c(k),'--']); 
grid on;
xlim([48,432]); ylim([-1.5 5.])
% legend(EF, 'EF1','EF2','EF3','EF4','EF5','EF6')
% legend('Orientation','horizontal')
% legend('Location','northoutside')
%legend('boxoff')
%title('FE com ruido SNR = 60 dB')
xlabel('\tau_n [amostras]')
ylabel('K_fE [Hz]')
title(tits(k));
end
% --------------------------------------------------
subplot(3,2,4); %ylim([-1.1 -0.9])
savefig(fig3,'salto_freq_KFE_tau_n.fig')
saveas(fig3,'salto_freq_KFE_tau_n.png')