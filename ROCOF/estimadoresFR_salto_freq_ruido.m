%testes para estimação de fase, frequência e ROCOF

clear all; close all; clc
%fixed parameters
SNR = 60; %Signal noise ratio[dB]
F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency
k_a = 0.0; % phase (angle) step height [degrees]
k_x = 0.0; % magnitude step height [relative]
k_f = 1.0; % frequency step height [Hz]
Fs = 4800; % sampling frequency [Hz]
NCycles = 6; % signal number of generated nominal cycles
T = NCycles/F0; % fundamental cycles duration
NSamples = floor(NCycles*Fs/F0); % total number of signal samples
phi_0 = 0; % initial angle phi_0 [degrees]
tau1 = 0.5;  % first step location in proportion of T
tau2 = 1.; % second step location in proportion of T; set tau2 = 1 if you dont want two steps
tau_n1 = floor(tau1*NSamples); %first step sample location
tau_n2 = floor(tau2*NSamples); %2nd step sample location
nbits = 32; % number of bits for the simulated signal generator
 
phistep = 5; %[degrees]
ncurves = 180/phistep;
phi_n = (70:phistep:(ncurves-1)*phistep) + phi_0;
Fref = (tau1*T*F1 + (T - tau1*T)*(F1 + k_f))/T; %one step only
ROCOF_ref = (Fref - F1)/T;

for j = 1:ncurves
    MC_iterations = 1;
    SNR_high = 1000; %sem ruido
    [F_est,f1e,f2e,kfe] = MC_estimation(MC_iterations,F0,F1,Fs,phi_n(j),NCycles,tau1,tau2,SNR_high,k_a, k_x,k_f,nbits);
    
    MC_iterations = 300;
    [Fraw,f1raw,f2raw,kfraw] = MC_estimation(MC_iterations,F0,F1,Fs,phi_n(j),NCycles,tau1,tau2,SNR,k_a, k_x,k_f,nbits);
    F_mean = mean(Fraw,2);
    F_mode = mode(Fraw,2);
    f1_mean(j,:) = mean(f1raw,2);
    f1_std(j,:) = std(f1raw,0,2);
    f2_mean(j,:) = mean(f2raw,2);
    f2_std(j,:) = std(f2raw,0,2);
    kf_mean(j,:) = mean(kfraw,2);
    kf_std(j,:) = std(f2raw,0,2);
    if MC_iterations>1
        FE_std(j,:) = std(Fraw,0,2);
    else
        FE_std(j,:) = 0;
    end
    
    % 1st pearson coeficiente de assimetria
    As(j,:) = (F_mean - F_mode)./FE_std(j,:)';
    
    FE(j,:) = F_est - Fref;
    FE_ruido(j,:) = F_mean - Fref;
        
    ROCOF_est = (F_est - Fref)./T;
    RFE(j,:) = ROCOF_est - ROCOF_ref; % [Hz/s]
end
FE
RFE

% figure
% plot(phi_n,(FE-mean(FE))/F1)
% xlabel('\phi_0 [graus]')
% ylabel('FE - \mu(FE) [Hz/Hz]')
% legend('EF1','EF2','EF3','EF4')
% 

c = ['k','b','c','r','m'];
fig1 = figure(1)
for k =1:size(FE,2)
plot(phi_n,FE(:,k),c(k)); hold on;
end
title('FE sem adicao de ruido')
xlabel('\phi_0 [graus]')
ylabel('FE [Hz]')
legend('EF1','EF2','EF3','EF4','EF5')

fig2 = figure(2)

for k = 1:size(FE,2)
EF(k) = plot(phi_n,FE_ruido(:,k),c(k)); hold on;
plot(phi_n,FE_ruido(:,k) + FE_std(:,k),[c(k),'--']);
plot(phi_n,FE_ruido(:,k) - FE_std(:,k),[c(k),'--']); 
end
title('FE com ruido SNR = 60 dB')
xlabel('\phi_0 [graus]')
ylabel('FE [Hz]')
legend(EF,'EF1','EF2','EF3','EF4','EF5','EF6')

ax1 = get(fig1,'CurrentAxes')
figure(2); ylim(ax1.YLim)

figure
FE_diff = FE - FE_ruido;
for k=1:size(FE,2)
plot(phi_n,FE_diff(:,k),c(k)); hold on;
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

fig5 = figure(5)

for k = 1:size(FE,2)
EF(k) = plot(phi_n,f1_mean(:,k),c(k)); hold on;
plot(phi_n,f1_mean(:,k) + f1_std(:,k),[c(k),'--']);
plot(phi_n,f1_mean(:,k) - f1_std(:,k),[c(k),'--']); 
end
title('F1 com ruido SNR = 60 dB')
xlabel('\phi_0 [graus]')
ylabel('f_1 [Hz]')
legend(EF,'EF1','EF2','EF3','EF4','EF5','EF6')

fig6 = figure(6)

for k = 1:size(FE,2)
EF(k) = plot(phi_n,kf_mean(:,k),c(k)); hold on;
plot(phi_n,kf_mean(:,k) + kf_std(:,k),[c(k),'--']);
plot(phi_n,kf_mean(:,k) - kf_std(:,k),[c(k),'--']); 
end
title('kf com ruido SNR = 60 dB')
xlabel('\phi_0 [graus]')
ylabel('kf [Hz]')
legend(EF,'EF1','EF2','EF3','EF4','EF5','EF6')

fig7 = figure(7)

for k = 1:size(FE,2)
EF(k) = plot(phi_n,f2_mean(:,k),c(k)); hold on;
plot(phi_n,f2_mean(:,k) + f2_std(:,k),[c(k),'--']);
plot(phi_n,f2_mean(:,k) - f2_std(:,k),[c(k),'--']); 
end
title('f_2 com ruido SNR = 60 dB')
xlabel('\phi_0 [graus]')
ylabel('f_2 [Hz]')
legend(EF,'EF1','EF2','EF3','EF4','EF5','EF6')


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