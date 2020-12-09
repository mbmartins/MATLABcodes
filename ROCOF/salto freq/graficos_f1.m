% influencia de tau
load('estimadoresFR_salto_freq_F1.mat');
close all;

c = ['k','b','c','r','m','g'];
tits = ["a)EF1","b)EF2","c)EF3","d)EF4","e)EF5","f)EF6"];

% ------------- FE com ruído de 60 dB -------------
fig1 = figure('Units','normalized','Position',[0 0 0.5 1])
for k = 1:6
subplot(3,2,k)
EF(k) = plot(F1_vec,FE(:,k),c(k)); hold on;
shade(F1_vec,FE(:,k) - FE_std(:,k),[c(k),'-.'],F1_vec,FE(:,k) + FE_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
%plot(tau_n,FE(:,k) - FE_std(:,k),[c(k),'--']); 
grid on;
%xlim([48,432]); 
ylim([-0.4 0.2])
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
savefig(fig1,'salto_freq_FE_F1.fig')
saveas(fig1,'salto_freq_FE_F1.png')

fig2 = figure('Units','normalized','Position',[0 0 0.5 1])
% ------------- FE1 com ruído de 60 dB -------------
for k = 1:6
subplot(3,2,k)
EF(k) = plot(F1_vec,FE1(:,k),c(k)); hold on;
shade(F1_vec,FE1(:,k) - f1_std(:,k),[c(k),'-.'],F1_vec,FE1(:,k) + f1_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
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
savefig(fig2,'salto_freq_F1_F1.fig')
saveas(fig2,'salto_freq_F1_F1.png')

fig3 = figure('Units','normalized','Position',[0 0 0.5 1]);
% ------------- KFE1 com ruído de 60 dB -------------
for k = 1:6
subplot(3,2,k)
EF(k) = plot(F1_vec,KFE(:,k),c(k)); hold on;
shade(F1_vec,KFE(:,k) - kf_std(:,k),[c(k),'-.'],F1_vec,KFE(:,k) + kf_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
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
savefig(fig3,'salto_freq_KFE_F1.fig')
saveas(fig3,'salto_freq_KFE_F1.png')
% 
% fig4 = figure('Units','normalized','Position',[0 0 0.5 0.5]);
% ------------- ALL com ruído de 60 dB -------------
for k = 1:6
    figALL(k) = figure('Units','normalized','Position',[0 0 0.5 0.13]);
    subplot(1,3,1)
    pleg(1) = plot(F1_vec,FE(:,k),c(k)+"-"); hold on;
    shade(F1_vec,FE(:,k) - FE_std(:,k),[c(k),'-.'],F1_vec,FE(:,k) + FE_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel('f_1 [Hz]'); ylabel('FE [Hz]')
    ylim([-0.3 0.1])
grid on;
    subplot(1,3,2)
    pleg(2) = plot(F1_vec,KFE(:,k),c(k)+"-"); hold on;
    shade(F1_vec,KFE(:,k) - kf_std(:,k),[c(k),'--'],F1_vec,KFE(:,k) + kf_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel('f_1 [Hz]'); ylabel('K_fE [Hz]')
    ylim([-1.1 0.1])
    grid on;
    subplot(1,3,3)
    pleg(3) = plot(F1_vec,FE1(:,k),c(k)+"-"); hold on;
    shade(F1_vec,FE1(:,k) - f1_std(:,k),[c(k),'--'],F1_vec,FE1(:,k) + f1_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel('f_1 [Hz]'); ylabel('FE_1 [Hz]')
    grid on;
    ylim([-0.3 0.6])
    %legend(pleg,'FE');
%plot(tau_n,FE(:,k) - FE_std(:,k),[c(k),'--']); 

%xlim([48,432]); 
%ylim([-0.4 0.1])
% legend(EF, 'EF1','EF2','EF3','EF4','EF5','EF6')
% legend('Orientation','horizontal')
% legend('Location','northoutside')
%legend('boxoff')
%title('FE com ruido SNR = 60 dB')
%title(tits(k));
savefig(figALL(k),"salto_freq_ALL_F1"+"_EF"+k)
saveas(figALL(k),"salto_freq_ALL_F1"+"_EF"+k+".png")
end
% --------------------------------------------------
