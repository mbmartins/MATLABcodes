%testes para estimação de fase, frequência e ROCOF

clear all; close all; clc
%fixed parameters
SNR = 60; %Signal noise ratio[dB]

%%% ajustar os valores de k para controlar os saltos
k_a = 10.0; % phase (angle) step height [degrees]
k_x = 0.0; % magnitude step height [relative]
k_f = 0.0; % frequency step height [Hz]
k_r = 0.0; % rocof step (not yet working)

F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency

Fs = 4800; % sampling frequency [Hz]
NCycles = 6; % signal number of generated nominal cycles
T = NCycles/F0; % fundamental cycles duration
NSamples = floor(NCycles*Fs/F0); % total number of signal samples
phi_0 = 0; % initial angle phi_0 [degrees]
tau1 = 0.5;  % first step location in proportion of T
tau2 = 1.; % second step location in proportion of T; set tau2 = 1 if you dont want two steps
tau_n1 = floor(tau1*NSamples); %first step sample location
tau_n2 = floor(tau2*NSamples); %2nd step sample location
nbits = 16; % number of bits for the simulated signal generator
 
phistep = 5; %[degrees]
ncurves = 180/phistep;
phi_n = (0:phistep:(ncurves-1)*phistep) + phi_0;
Fref = (tau1*T*F1 + (T - tau1*T)*(F1 + k_f))/T; %one step only
ROCOF_ref = (Fref - F1)/T;

for j = 1:ncurves
    MC_iterations = 1;
    SNR_high = 1000; %sem ruido
    [F_est,f1e,f2e,kfe] = MC_estimation(MC_iterations,F0,F1,Fs,phi_n(j),NCycles,tau1,tau2,SNR_high,k_a, k_x,k_f,nbits);
    
    MC_iterations = 300;
    [Fraw,f1raw,f2raw,kfraw] = MC_estimation(MC_iterations,F0,F1,Fs,phi_n(j),NCycles,tau1,tau2,SNR,k_a, k_x,k_f,nbits);
    F_mean = mean(Fraw,2);
    %F_mode = mode(Fraw,2);
    
    f1_mean(j,:) = mean(f1raw,2);
    f1_std(j,:) = std(f1raw,0,2);
    
    f2_mean(j,:) = mean(f2raw,2);
    f2_std(j,:) = std(f2raw,0,2);
    kf_mean(j,:) = mean(kfraw,2);
    kf_std(j,:) = std(kfraw,0,2);
    if MC_iterations>1
        FE_std(j,:) = std(Fraw,0,2);
    else
        FE_std(j,:) = 0;
    end
    
    % 1st pearson coeficiente de assimetria
%    As(j,:) = (F_mean - F_mode)./FE_std(j,:)';
    
    FE(j,:) = F_est - Fref;
    FE_ruido(j,:) = F_mean - Fref;
   
    FE1(j,:) = f1e - F1;
    FE1_ruido(j,:) = f1_mean(j,:) - F1;
    
    KE(j,:) = kfe - k_f;
    KE_ruido(j,:) = kf_mean(j,:) - k_f;
    
    ROCOF_est = (F_est - Fref)./T;
    RFE(j,:) = ROCOF_est - ROCOF_ref; % [Hz/s]
end

% figure
% plot(phi_n,(FE-mean(FE))/F1)
% xlabel('\phi_0 [graus]')
% ylabel('FE - \mu(FE) [Hz/Hz]')
% legend('EF1','EF2','EF3','EF4')
% 

c = ['k','b','c','r','m','g'];
fig1 = figure(1)
%subplot(1,3,1)
for k =1:6
plot(phi_n,FE(:,k),c(k)); hold on;
end
title('FE sem adicao de ruido')
xlabel('\phi_0 [graus]')
ylabel('FE [Hz]')
legend('EF1','EF2','EF3','EF4','EF5','EF6')
% subplot(1,3,2)
% for k =4:6
% plot(phi_n,FE(:,k),c(k)); hold on;
% end
% xlabel('\phi_0 [graus]')
% ylabel('FE [Hz]')
% legend('EF4','EF5','EF6')
% subplot(1,3,3)
% plot(phi_n,FE(:,2),c(2)); hold on;
% plot(phi_n,FE(:,4),c(4));
% xlabel('\phi_0 [graus]')
% ylabel('FE [Hz]')
% legend('EF2','EF4')

fig2 = figure(2)

%subplot(1,3,1)
for k = 1:6
EF(k) = plot(phi_n,FE_ruido(:,k),c(k)); hold on;
plot(phi_n,FE_ruido(:,k) + FE_std(:,k),[c(k),'--']);
plot(phi_n,FE_ruido(:,k) - FE_std(:,k),[c(k),'--']); 
end
legend(EF, 'EF1','EF2','EF3','EF4','EF5','EF6')
%title('FE com ruido SNR = 60 dB')
xlabel('\phi_0 [graus]')
ylabel('FE [Hz]')
% subplot(1,3,2)
% for k = 4:6
% EFP(k-3) = plot(phi_n,FE_ruido(:,k),c(k)); hold on;
% plot(phi_n,FE_ruido(:,k) + FE_std(:,k),[c(k),'--']);
% plot(phi_n,FE_ruido(:,k) - FE_std(:,k),[c(k),'--']); 
% end
% xlabel('\phi_0 [graus]')
% ylabel('FE [Hz]')
% legend(EFP,'EF4','EF5','EF6')
% subplot(1,3,3)
% k = 2;
% EFp3(1) = plot(phi_n,FE_ruido(:,k),c(k)); hold on;
% plot(phi_n,FE_ruido(:,k) + FE_std(:,k),[c(k),'--']);
% plot(phi_n,FE_ruido(:,k) - FE_std(:,k),[c(k),'--']); 
% k = 4;
% EFp3(2) = plot(phi_n,FE_ruido(:,k),c(k)); hold on;
% plot(phi_n,FE_ruido(:,k) + FE_std(:,k),[c(k),'--']);
% plot(phi_n,FE_ruido(:,k) - FE_std(:,k),[c(k),'--']);
% xlabel('\phi_0 [graus]')
% ylabel('FE [Hz]')
% legend(EFp3,'EF2','EF4')


figure(3); hold off;
%s = ["kx",'b','c','r','m','g'];
s = c;
FE_diff = FE - FE_ruido;
for k=1:size(FE,2)
plot(phi_n,FE_diff(:,k),s(k)); hold on;
end
xlabel('\phi_0 [graus]')
ylabel('FE - FE_{ruido} [Hz]')
legend('EF1','EF2','EF3','EF4','EF5','EF6')
% 
% figure
% k=4
% plot(phi_n,FE_diff(:,k),c(k)); hold on;
% xlabel('\phi_0 [graus]')
% ylabel('FE - FE_{ruido} [Hz]')
% legend('EF4')

% ----- Grafico da estimacao de f1 ------%
fig5 = figure(5)
subplot(1,3,1)
for k = 1:6
plotf1(k) = plot(phi_n,FE1(:,k),c(k)); hold on;
%plot(phi_n,FE1(:,k) + f1_std(:,k),[c(k),'--']);
%plot(phi_n,FE1(:,k) - f1_std(:,k),[c(k),'--']); 
end
ylim([-0.1 0.1]);
title('a)');% SNR = 60 dB')
xlabel('\phi_0 [graus]')
ylabel('FE_1 [Hz]')
%legend(plotf1,'EF1','EF2','EF3','EF4','EF5','EF6')
subplot(1,3,2)
for k = 1:6
  plotf1P(k) = plot(phi_n,FE1_ruido(:,k),c(k)); hold on;
 plot(phi_n,FE1_ruido(:,k) + f1_std(:,k),[c(k),'--']);
 plot(phi_n,FE1_ruido(:,k) - f1_std(:,k),[c(k),'--']); 
end
ylim([-0.1 0.1]);
title('b)')
 xlabel('\phi_0 [graus]')
 ylabel('FE_1 [Hz]')
%legend(plotf1,'EF1','EF2','EF3','EF4','EF5','EF6')
% 
subplot(1,3,3)
FE1diff = FE1 - FE1_ruido;
for k = 1:6
  plotf1P(k) = plot(phi_n,FE1diff(:,k),c(k)); hold on;
end
ylim([-0.1 0.1]);
title('c)')
 xlabel('\phi_0 [graus]')
 ylabel('FE_1 - FE_1_{ruido} [Hz]')
legend(plotf1P,'EF1','EF2','EF3','EF4','EF5','EF6')




%  k = 2;hold off;
%  EFp3(1) = plot(phi_n,FE1(:,k),c(k)); hold on;
%  plot(phi_n,FE1(:,k) + f1_std(:,k),[c(k),'--']);
%  plot(phi_n,FE1(:,k) - f1_std(:,k),[c(k),'--']); 
%  k = 5;
%  EFp3(2) = plot(phi_n,FE1(:,k),c(k)); hold on;
% % plot(phi_n,FE1(:,k) + f1_std(:,k),[c(k),'--']);
% % plot(phi_n,FE1(:,k) - f1_std(:,k),[c(k),'--']);
% % xlabel('\phi_0 [graus]')
% % ylabel('FE_1 [Hz]')
% % legend(EFp3,'EF2','EF5')


%----- Grafico da estimacao de kf ----- %
fig6 = figure(6)
subplot(1,3,1)
for k = 1:6
plotk(k) = plot(phi_n,KE(:,k),c(k)); hold on;
%plot(phi_n,FE1(:,k) + f1_std(:,k),[c(k),'--']);
%plot(phi_n,FE1(:,k) - f1_std(:,k),[c(k),'--']); 
end
ylim([-0.15 0.15]);
title('a)');% SNR = 60 dB')
xlabel('\phi_0 [graus]')
ylabel('K_fE [Hz]')
%legend(plotf1,'EF1','EF2','EF3','EF4','EF5','EF6')
subplot(1,3,2)
for k = 1:6
  plotkP(k) = plot(phi_n,KE_ruido(:,k),c(k)); hold on;
 plot(phi_n,KE_ruido(:,k) + kf_std(:,k),[c(k),'--']);
 plot(phi_n,KE_ruido(:,k) - kf_std(:,k),[c(k),'--']); 
end
 title('b)')
 ylim([-0.15 0.15]);
 xlabel('\phi_0 [graus]')
 ylabel('K_fE [Hz]')
legend(plotkP,'EF1','EF2','EF3','EF4','EF5','EF6')
% 
subplot(1,3,3)
KEdiff = KE - KE_ruido;
for k = 1:6
  plotKEP(k) = plot(phi_n,KEdiff(:,k),c(k)); hold on;
end
ylim([-0.15 0.15]);
 title('c)')
 xlabel('\phi_0 [graus]')
 ylabel('K_fE - K_fE_{ruido} [Hz]')
legend(plotKEP,'EF1','EF2','EF3','EF4','EF5','EF6')


% 
% for k = 1:6
% EF(k) = plot(phi_n,kf_mean(:,k),c(k)); hold on;
% plot(phi_n,kf_mean(:,k) + kf_std(:,k),[c(k),'--']);
% plot(phi_n,kf_mean(:,k) - kf_std(:,k),[c(k),'--']); 
% end
% title('kf com ruido SNR = 60 dB')
% xlabel('\phi_0 [graus]')
% ylabel('kf [Hz]')
% legend(EF,'EF1','EF2','EF3','EF4','EF5','EF6')
% 
% fig7 = figure(7)
% 
% for k = 1:size(FE,2)
% EF(k) = plot(phi_n,f2_mean(:,k),c(k)); hold on;
% plot(phi_n,f2_mean(:,k) + f2_std(:,k),[c(k),'--']);
% plot(phi_n,f2_mean(:,k) - f2_std(:,k),[c(k),'--']); 
% end
% title('f_2 com ruido SNR = 60 dB')
% xlabel('\phi_0 [graus]')
% ylabel('f_2 [Hz]')
% legend(EF,'EF1','EF2','EF3','EF4','EF5','EF6')


% figure
% hist(Fraw(2,:))
% title('FE EF2')
% 
% figure
% plot(phi_n,phi_0_est)
% xlabel('\phi_0 [graus]')
% ylabel('phi_0_est [graus]')
% 
% figure
% fref = Fref*ones(size(f_u1));
% plot(f_u1); hold on; plot(f_u2); plot(f_u3); plot(f_u4);plot(fref,'k')
% legend('fu1','fu2','fu3','fu4')