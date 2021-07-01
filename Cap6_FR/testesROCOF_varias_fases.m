%testes para estimação de fase, frequência e ROCOF
clear all; close all; clc
SNR = 600;
%angle phi_0 [degrees]
phi_0 = 0;
%fixed parameters
F0 = 60;
F1 = 60;
tau1 = 0.7;  % in proportion of the time window
tau2 = 2; % in proportion of the time window
KaS = 10.0; %[degrees]
KxS = 0.0; % [relative magnitude step]
KfS = 0.; %[Hz] %size of frequency step
Fs =4800;
NCycles = 6;
T = NCycles/F0;
NSamples = floor(NCycles*Fs/F0);
tau_n = floor(tau1*NSamples)

phistep = 5; %[degrees]
ncurves = 360/phistep;
phi_n = (0:phistep:(ncurves-1)*phistep) + phi_0;
n = 1:NSamples;
br = 0.05; brn = br*NSamples; %indexes of samples to ignore
br_mask = floor((br)*NSamples):floor((1-br)*NSamples);

for j = 1:ncurves
    Signal(j,:) = SigGEN2(F0,F1,Fs,phi_n(j),NCycles,tau1,tau2,SNR,KaS, KxS,KfS);
    z=hilbert(Signal(j,:)');
    Psi_i(:,j) = unwrap(angle(z)); %[rad]
    %Psi_i = phase(z);
    dfh(:,j)=gradient(Psi_i(:,j))*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    mz = abs(z);
    df(:,j) = dfh(:,j).*mz./median(mz); %compensation of magnitude interference
    fest(j) = median(dfh(:,j));
    fest_comp(j) = median(df(br_mask,j));
    ROCOF_ih(:,j) = gradient(dfh(:,j))/T;
    ROCOF_i(:,j) = gradient(df(:,j))/T;
    ROCOF_est(j) = median(ROCOF_ih(:,j));
    ROCOF_est_comp(j) = median(ROCOF_i(:,j))
    f1est(j) = median(df(1:tau_n,j));
    f2est(j) = median(df(tau_n+1:NSamples,j));   
    fest_compDF2(j) = (f1est(j) + f2est(j))/2;
    ROCOF_estDF2(j) = (f1est(j) - f2est(j))/T;
end


%% Figures for Psi[n] signals
figure(1);
for j = 1:ncurves
    plot(n(225:255),Psi_i(225:255,j)*180/pi);hold on;
end
xlabel('Amostras'); 
ylabel('\Psi[n] [°]')

% figure(2)
%  surf(gradient(Psi_i(225:255,:)))
 
 figure(2)
 plot(phi_n,Psi_i(241,:)*180/pi)
xlabel('\phi_0'); 
ylabel('\Psi[\tau_n] [°]')

figure(3)
plot(phi_n,df(tau_n+2,:),'--'); hold on;
plot(phi_n,df(tau_n+1,:),'--')
plot(phi_n,df(tau_n,:))
plot(phi_n,df(tau_n-1,:),'--')
plot(phi_n,df(tau_n-2,:),'--')
xlabel('\phi_0 [°]'); 
ylabel('f_i[\tau_n] [Hz]')
legend('f_i[\tau_n+2]','f_i[\tau_n+1]','f_i[\tau_n]','f_i[\tau_n-1]','f_i[\tau_n-2]' )

figure(4)
plot(phi_n,ROCOF_i(tau_n+2,:),'--'); hold on;
plot(phi_n,ROCOF_i(tau_n+1,:),'--')
plot(phi_n,ROCOF_i(tau_n,:))
plot(phi_n,ROCOF_i(tau_n-1,:),'--')
plot(phi_n,ROCOF_i(tau_n-2,:),'--')
xlabel('\phi_0 [°]'); 
ylabel('ROCOF_i[n] [Hz/s]')
legend('ROCOF_i[\tau_n+2]','ROCOF_i[\tau_n+1]','ROCOF_i[\tau_n]','ROCOF_i[\tau_n-1]','ROCOF_i[\tau_n-2]' )

f1vec = F1*ones(size(phi_n));
figure;
plot(phi_n,fest,'.-',phi_n,fest_comp,'.-',phi_n,fest_compDF2,phi_n,f1vec,'-')
title('Frequência estimada em função de \phi_0 para salto de fase de 10 graus')
xlabel('\phi_0 [°]'); 
ylabel('f [Hz]')
legend('f sem compensação','f com compensação','f por DF2','f_1')

figure;
%ROCOF_vec = [0.5]'*ones(size(phi_n));
plot(phi_n,ROCOF_est,'.-',phi_n,ROCOF_est_comp,'.-',phi_n,ROCOF_estDF2) %,phi_n,ROCOF_vec,'k',phi_n,-ROCOF_vec,'k')
title('ROCOF de rede estimada em função de \phi_0 para salto de fase de 10 graus')
xlabel('\phi_0 [°]'); 
ylabel('ROCOF [Hz/s]')
legend('ROCOF sem compensação','ROCOF com compensação','ROCOF DF2')

% 
% figure;
% % o erro absoluto
% title('Erros de frequência')
% FE = (fest_comp - F1); 
% plot(phi_n,FE,'.-');
% xlabel('\phi_0 [°]'); 
% ylabel('FE [Hz]')
% 
% figure;histfit(FE)
% FE_mean = mean(FE)
% FE_std = std(FE)