% influencia de phi_0
load('salto_freq_phi0.mat');
close all;

c = ['k','b','c','r','m','g'];
tits = ["a)EF1","b)EF2","c)EF3","d)EF4","e)EF5","f)EF6"];

% ------------- FE com ruído de 60 dB -------------
fig1 = figure(1)

for k = 1:6
subplot(3,2,k)
xlim([0,180]); ylim([-0.06 0.05])
EF(k) = plot(phi_n,FE_ruido(:,k),c(k)); hold on;
shade(phi_n,FE_ruido(:,k) - FE_std(:,k),[c(k),'-.'],phi_n,FE_ruido(:,k) + FE_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
%plot(phi_n,FE_ruido(:,k) - FE_std(:,k),[c(k),'--']); 
grid on;
xlim([0,175]); ylim([-0.06 0.06])
% legend(EF, 'EF1','EF2','EF3','EF4','EF5','EF6')
% legend('Orientation','horizontal')
% legend('Location','northoutside')
%legend('boxoff')
%title('FE com ruido SNR = 60 dB')
xlabel('\phi_0 [graus]')
ylabel('FE [Hz]')
title(tits(k));
end
% --------------------------------------------------
savefig(fig1,'salto_freq_phi0.fig')

fig2 = figure(2)
% ------------- FE1 com ruído de 60 dB -------------
for k = 1:6
subplot(3,2,k)
EF(k) = plot(phi_n,FE1_ruido(:,k),c(k)); hold on;
shade(phi_n,FE1_ruido(:,k) - f1_std(:,k),[c(k),'-.'],phi_n,FE1_ruido(:,k) + f1_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
%plot(phi_n,FE_ruido(:,k) - FE_std(:,k),[c(k),'--']); 
grid on;
xlim([0,175]); ylim([-0.1 0.1])
% legend(EF, 'EF1','EF2','EF3','EF4','EF5','EF6')
% legend('Orientation','horizontal')
% legend('Location','northoutside')
%legend('boxoff')
%title('FE com ruido SNR = 60 dB')
xlabel('\phi_0 [graus]')
ylabel('FE1 [Hz]')
title(tits(k));
end
% --------------------------------------------------
subplot(3,2,4); ylim([0.4 0.5])
savefig(fig2,'salto_freq_F1_phi0.fig')

fig3 = figure(3);
% ------------- KFE1 com ruído de 60 dB -------------
for k = 1:6
subplot(3,2,k)
EF(k) = plot(phi_n,KE_ruido(:,k),c(k)); hold on;
shade(phi_n,KE_ruido(:,k) - kf_std(:,k),[c(k),'-.'],phi_n,KE_ruido(:,k) + kf_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
%plot(phi_n,FE_ruido(:,k) - FE_std(:,k),[c(k),'--']); 
grid on;
xlim([0,175]); ylim([-0.2 0.2])
% legend(EF, 'EF1','EF2','EF3','EF4','EF5','EF6')
% legend('Orientation','horizontal')
% legend('Location','northoutside')
%legend('boxoff')
%title('FE com ruido SNR = 60 dB')
xlabel('\phi_0 [graus]')
ylabel('K_fE [Hz]')
title(tits(k));
end
% --------------------------------------------------
subplot(3,2,4); ylim([-1.1 -0.9])
savefig(fig3,'salto_freq_KFE_phi0.fig')
