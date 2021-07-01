%testes para estimação de fase, frequência e ROCOF
clear all; close all; clc
SNR = 600;
%angle phi_0 [degrees]
phi_0 = 175;
%fixed parameters
F0 = 60.0; %nominal frequency
F1 = F0; %fundamental frequency
tau1 = 0.5;  % in proportion of the time window
tau2 = 2; % in proportion of the time window
KaS = 10.0; %[degrees]
KxS = 0.0; % [relative magnitude step]
KfS = 0.; %[Hz] %size of frequency step
Fs =4800;
NCycles = 6;
T = NCycles/F0;
NSamples = floor(NCycles*Fs/F0);
tau_n = floor(tau1*NSamples)

Signal = SigGEN2(F0,F1,Fs,phi_0,NCycles,tau1,tau2,SNR,KaS, KxS,KfS);
z=hilbert(Signal');
Psi_i = unwrap(angle(z)); %[rad]
%Psi_i = phase(z);

dfh=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
mz = abs(z);
df = dfh.*mz./median(mz);
ROCOF_ih = gradient(dfh)/T;  %% verificar se precisa ponderar por T
ROCOF_i = gradient(df)/T;    %% na verdade nao iremos usar estes sinais

F_est[1] = F1;
ROCOF_est = median(ROCOF_i);

% Estimation and analysis
n = 1:NSamples;

%EF1
fu_est = 


br = 0.05; brn = br*NSamples; %indexes of samples to ignore
br_mask = floor((br)*NSamples):floor((1-br)*NSamples);

fest = median(df(brn:end-brn)); %frequency estimate
%fest = median(df(50:100)); %frequency estimate

%plot(Signal);
%acertar os eixos x

F1v(1:NSamples) = F1;
FE = 100*(fest - F1)/F1
dferr = df - F1v;
% figure; plot(dferr); title('fi error')

%%PATV filter for Psi
lambda = 1; d = 1; Nit = 30;
[x, p, cost, ruido, v] = patv_MM(Psi_i, d, lambda, Nit);
Psi_PATV = p + x;
df1 = gradient(Psi_PATV)*Fs/(2*pi); %[Hz]


%% Figures for Psi[n] signals
figure
plot(Psi_PATV);
title('PATV filtered \Psi[n]')
xlabel('Amostras'); ylabel('\Psi[n] [rad/s]')

figure
plot(dfh); hold on; plot(df); %plot(df1)
xlabel('Amostras'); ylabel('f_i [Hz]')
legend('f_i sem compensação','f_i com compensação')%,'f_i 1')

figure
plot(ROCOF_ih); hold on; plot(ROCOF_i)
xlabel('Amostras'); ylabel('ROCOF_i [Hz/s]')
legend('ROCOF_i sem compensação','ROCOF_i com compensação')


figure(2);
%subplot(3,1,1)
%plot(n(1:tau_n),Psi_i(1:tau_n)*180/pi,'b'); hold on;
%xlabel('Amostras'); 
%plot(n(tau_n:end),Psi_i(tau_n:end)*180/pi,'r');
%xlabel('Amostras'); 
%ylabel('\Psi[n] [°]')
%legend('\Psi_1[n]','\Psi_2[n]')
subplot(3,1,1)
plot(n(tau_n-10:tau_n),Psi_i(tau_n-10:tau_n)*180/pi,'b'); hold on;
%xlabel('Amostras');
plot(n(tau_n:tau_n+10),Psi_i(tau_n:tau_n+10)*180/pi,'r');
%xlabel('Amostras'); 
ylabel('\Psi[n] [°]')
legend('\Psi_1[n]','\Psi_2[n]')

%hold on; plot(Psi_PATV); legend('\Psi','\Psi_{PATV}')
subplot(3,1,2)
deltaPsi = gradient(Psi_i)*180/pi;
%plot(v);hold on; legend('v ~ f_i');
plot(n(tau_n-10:tau_n+10),deltaPsi(tau_n-10:tau_n+10)); 
ylabel('\Delta\Psi[n] [°]')
xlabel('Amostras');

subplot(3,1,3)
%plot(v);hold on; legend('v ~ f_i');
plot(n(tau_n-10:tau_n+10),df(tau_n-10:tau_n+10)); 
ylabel('f_i [Hz]')
xlabel('Amostras');

% 
% hold on; plot(gradient(Psi_PATV));
% legend('\Delta \Psi','\Delta \Psi_{PATV}')

%% PATV filter for fi
lambda = 3; d = 1; Nit = 30;
[x, p, cost, u, v] = patv_MM(df, d, lambda, Nit);
fi_PATV = p+x;
ROCOF_PATV = u; %gradient(p+x);
figure; plot(fi_PATV); title('PATV-filtered f_i')
hold on; plot(df); %plot(ROCOF_PATV); 
legend('f_i-PATV','f_i'); %,'ROCOF_i')
figure; plot(ROCOF_i); 
hold on; plot(ROCOF_PATV); title('estimated ROCOF_i')
plot(v); legend('ROCOF_i','ROCOF_i-PATV','extracted noise')

% o PATV pode ser especialmente útil na discriminação do salto de fase
% combinado com o salto de frequência. Nos casos individuais não é de tanta
% serventia.

%podemos desenvolver um método de baixo custo computacional, mas com erros
%mais altos, usando técnicas mais simples; e outro que use mais iterações,
%PATV, etc, para atingir maior exatidão e robustez.

%%% identificar os limites de f1 e f2 a partir da detecção de tau
% ignora-se 5% amostras prox as bordas
f1 = (df((1:tau_n-brn)));
f2 = (df(tau_n-brn:end));

% figure; histogram(f1,'BinLimits',[55,75])
% figure; histogram(f2,'BinLimits',[55,75])

f1m = median(f1);
f2m = median(f2);
ROCOF = f2m - f1m
rcum = sum(ROCOF_i(tau_n-2:tau_n+2))
% esta forma de identificar o ROCOF é tão 'preditiva' quanto a que usa
% modelo de Taylor, com a vantagem de ser direta.
% -> pode ser ponderada pela duração de tau, como o fasor intermediário
% -> desenvolver o sinal de detecção baseado em ROCOF para tau, como no paper
% do AMPS
% -> estabelecer critérios de detecção para sistema de proteção baseado no
% ERAC, do procedimentos de rede da ANEEL
% -> haverá uma faixa 'indeterminada', que depende da incerteza de medição

