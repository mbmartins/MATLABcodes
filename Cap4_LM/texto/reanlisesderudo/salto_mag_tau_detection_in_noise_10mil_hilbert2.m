clear all; close all; clc
fs=5000;% frequencia de amostragem 
f0=60; % frequencia nominal em Hz
A=10; % magnitude nominal 
fasei=90; % fase inicial em graus
phi=fasei*pi/180; % fase inicial em rad 
%N=4096; p=0.08157; 
N=500; 
pertau=.5;  % localização do salto em termos do percentual do tamanho da janela 
%tau=round(N*pertau); % índice do instante (arredondado) de ocorrência do salto
tau=292;  % tau ajustado para o salto de fase coincidir com um maximo local do seno
%tau=272; % tau ajustado para o salto de fase coincidir com um minimo local do seno
n=(0:(N-1))'; % base de tempo unitário
fase=2*pi*f0.*n/fs+phi; % fase limpa
SNRdB=51; % SRN alvo em dB
mag=[A*ones(tau-1,1); 1.1*A*ones(N-tau+1,1)]; 
x=mag.*cos(fase); % sinal sem ruído

vx=var(x); % variância do ruído
cn=sqrt(vx/10^(SNRdB/10))*randn(length(x),10000); % componentes de ruído
for k=1:10000,
xn=x+cn(:,k); % sinal com ruído
xa=hilbert(xn); 
% detetor 1
phase_est=phase(xa);
fi_anomaly=gradient(phase_est);
fi_anomaly=fi_anomaly-median(fi_anomaly(10:end-10));
fi_anomaly=fi_anomaly(10:end-10);
[v1,tau_est1]=max(abs(fi_anomaly));
limiar1=7*median(abs(fi_anomaly));

% detetor 2
mag_est=abs(xa);
mag_anomaly=gradient(mag_est);
mag_anomaly=mag_anomaly-median(mag_anomaly(10:end-10));
mag_anomaly=mag_anomaly(10:end-10);
limiar2=7*median(abs(mag_anomaly));
[v2,tau_est2]=max(abs(mag_anomaly));
if v1>limiar1,
    tau_est_v(k)=tau_est1(1)+9;
elseif v2>limiar2,
    tau_est_v(k)=tau_est2(1)+9;
else
    tau_est_v(k)=NaN;
end
end
hist (tau_est_v,tau-100:tau+100)

