% influencia de tau
load('estimadoresFR_salto_freq_kf.mat');
close all;

c = ['k','b','c','r','m','g'];
tits = ["a)EF1","b)EF2","c)EF3","d)EF4","e)EF5","f)EF6"];

% ------------- FE com ruído de 60 dB -------------
fig1 = figure('Units','normalized','Position',[0 0 0.5 1])
for k = 1:6
subplot(3,2,k)
EF(k) = plot(kf_vec,FE(:,k),c(k)); hold on;
shade(kf_vec,FE(:,k) - FE_std(:,k),[c(k),'-.'],kf_vec,FE(:,k) + FE_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
%plot(tau_n,FE(:,k) - FE_std(:,k),[c(k),'--']); 
grid on;
%xlim([48,432]); 
ylim([-0.3 0.1])
% legend(EF, 'EF1','EF2','EF3','EF4','EF5','EF6')
% legend('Orientation','horizontal')
% legend('Location','northoutside')
%legend('boxoff')
%title('FE com ruido SNR = 60 dB')
xlabel('f_1 [Hz]')
ylabel('FE [Hz]')
title(tits(k));
end
% --------------------------------------------------
savefig(fig1,'salto_freq_FE_kf.fig')
saveas(fig1,'salto_freq_FE_kf.png')

fig2 = figure('Units','normalized','Position',[0 0 0.5 1])
% ------------- FE1 com ruído de 60 dB -------------
for k = 1:6
subplot(3,2,k)
EF(k) = plot(kf_vec,FE1(:,k),c(k)); hold on;
shade(kf_vec,FE1(:,k) - f1_std(:,k),[c(k),'-.'],kf_vec,FE1(:,k) + f1_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
%plot(tau_n,FE(:,k) - FE_std(:,k),[c(k),'--']); 
grid on;
%xlim([48,432]); 
ylim([-0.3 0.6])
% legend(EF, 'EF1','EF2','EF3','EF4','EF5','EF6')
% legend('Orientation','horizontal')
% legend('Location','northoutside')
%legend('boxoff')
%title('FE com ruido SNR = 60 dB')
xlabel('f_1 [Hz]')
ylabel('FE1 [Hz]')
title(tits(k));
end
% --------------------------------------------------
%subplot(3,2,4); ylim([0.4 0.5])
savefig(fig2,'salto_freq_F1_kf.fig')
saveas(fig2,'salto_freq_F1_kf.png')

fig3 = figure('Units','normalized','Position',[0 0 0.5 1]);
% ------------- KFE1 com ruído de 60 dB -------------
for k = 1:6
subplot(3,2,k)
EF(k) = plot(kf_vec,KFE(:,k),c(k)); hold on;
shade(kf_vec,KFE(:,k) - kf_std(:,k),[c(k),'-.'],kf_vec,KFE(:,k) + kf_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
%plot(tau_n,FE(:,k) - FE_std(:,k),[c(k),'--']); 
grid on;
%xlim([48,432]); 
ylim([-0.4 0.1])
% legend(EF, 'EF1','EF2','EF3','EF4','EF5','EF6')
% legend('Orientation','horizontal')
% legend('Location','northoutside')
%legend('boxoff')
%title('FE com ruido SNR = 60 dB')
xlabel('f_1 [Hz]')
ylabel('K_fE [Hz]')
title(tits(k));
end
% --------------------------------------------------
subplot(3,2,4); ylim([-1.25 -0.75])
savefig(fig3,'salto_freq_KFE_kf.fig')
saveas(fig3,'salto_freq_KFE_kf.png')

fig4 = figure('Units','normalized','Position',[0 0 0.5 1]);
% ------------- ALL com ruído de 60 dB -------------
for k = 1:6
subplot(3,2,k)
pleg(1) = plot(kf_vec,FE(:,k),c(k)+"-"); hold on;
pleg(2) = plot(kf_vec,KFE(:,k),c(k)+"--"); hold on;
pleg(3) = plot(kf_vec,FE1(:,k),c(k)+".-"); hold on;
shade(kf_vec,FE(:,k) - FE_std(:,k),[c(k),'-.'],kf_vec,FE(:,k) + FE_std(:,k),[c(k),'-'], 'FillType', [1,2; 2,1]);
shade(kf_vec,KFE(:,k) - kf_std(:,k),[c(k),'-.'],kf_vec,KFE(:,k) + kf_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
shade(kf_vec,KFE(:,k) - f1_std(:,k),[c(k),'-.'],kf_vec,KFE(:,k) + f1_std(:,k),[c(k),'.-'], 'FillType', [1,2; 2,1]);
legend(pleg,'FE','K_fE','FE_1');
%plot(tau_n,FE(:,k) - FE_std(:,k),[c(k),'--']); 
grid on;
%xlim([48,432]); 
%ylim([-0.4 0.1])
% legend(EF, 'EF1','EF2','EF3','EF4','EF5','EF6')
% legend('Orientation','horizontal')
% legend('Location','northoutside')
%legend('boxoff')
%title('FE com ruido SNR = 60 dB')
xlabel('f_1 [Hz]')
ylabel('Error [Hz]')
title(tits(k));
end
% --------------------------------------------------
savefig(fig4,'salto_freq_ALL_kf.fig')
saveas(fig4,'salto_freq_ALL_kf.png')