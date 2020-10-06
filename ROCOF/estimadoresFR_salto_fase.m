%testes para estimação de fase, frequência e ROCOF

clear all; close all; clc
%fixed parameters
F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency
KaS = 10.0; %[degrees]
KxS = 0.0; % [relative magnitude step]
KfS = 0.0; %[Hz] %size of frequency step
Fs =4800;
NCycles = 6;
T = NCycles/F0;
phi_0 = 0; %angle phi_0 [degrees]
tau1 = 0.5;  % in proportion of T
tau2 = 1.; % in proportion of T; set tau2 = 1 if you dont want two steps
NSamples = floor(NCycles*Fs/F0);
tau_n1 = floor(tau1*NSamples);
tau_n2 = floor(tau2*NSamples);
nbits = 16;

phistep = 5; %[degrees]
ncurves = 180/phistep;
phi_n = (0:phistep:(ncurves-1)*phistep) + phi_0;
Fref = (tau1*T*F1 + (T - tau1*T)*(F1 + KfS))/T; %one step only
ROCOF_ref = (Fref - F1)/T;

for j = 1:ncurves
    MC_iterations = 1;
    SNR = 1000; %sem ruido
    F_est = MC_estimation(MC_iterations,F0,F1,Fs,phi_n(j),NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits);
    
    MC_iterations = 300;
    SNR = 60; % ruido de 60 dB
    Fraw = MC_estimation(MC_iterations,F0,F1,Fs,phi_n(j),NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits);
    F_mean = mean(Fraw,2);
    if MC_iterations>1
        FE_std(j,:) = std(Fraw,0,2);
    else
        FE_std(j,:) = 0;
    end
    
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

c = ['k','b','g','r'];
fig1 = figure(1)
for k =1:4
plot(phi_n,FE(:,k),c(k)); hold on;
end
xlabel('\phi_0 [graus]')
ylabel('FE [Hz]')
legend('EF1','EF2','EF3','EF4')

fig2 = figure(2)
for k = 1:4
EF(k) = plot(phi_n,FE_ruido(:,k),c(k)); hold on;
plot(phi_n,FE_ruido(:,k) + FE_std(:,k),[c(k),'--']);
plot(phi_n,FE_ruido(:,k) - FE_std(:,k),[c(k),'--']); 
end
xlabel('\phi_0 [graus]')
ylabel('FE [Hz]')
legend(EF,'EF1','EF2','EF3','EF4')

ax2 = get(fig2,'CurrentAxes')
figure(1); ylim(ax2.YLim)

figure;
subplot(1,2,1)
show(fig1)

figure
FE_diff = FE - FE_ruido;
for k=1:4
plot(phi_n,FE_diff(:,k),c(k)); hold on;
end
xlabel('\phi_0 [graus]')
ylabel('FE - FE_{ruido} [Hz]')
legend('EF1','EF2','EF3','EF4')

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