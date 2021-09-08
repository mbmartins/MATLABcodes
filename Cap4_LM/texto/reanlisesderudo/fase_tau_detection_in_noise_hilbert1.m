clear all; close all; clc
fs=5000;% frequencia de amostragem 
f0=60; % frequencia nominal em Hz
A=10; % magnitude nominal 
fasei=90; % fase inicial em graus
phi=fasei*pi/180; % fase inicial em rad 
N=500; 
pertau=.5;  % localização do salto em termos do percentual do tamanho da janela 
%tau=round(N*pertau); % índice do instante (arredondado) de ocorrência do salto
tau=312;  % tau ajustado para o salto de fase coincidir com um maximo local do seno
n=(0:(N-1))'; % base de tempo unitário
fase=2*pi*f0.*n/fs+phi; % fase limpa
fase(tau:end)=fase(tau:end)+10*pi/180; % salto de 10 graus na fase em t=tau (conforme eq (2))
%fase(tau:end)=fase(tau:end)+10*pi/180-(fase(tau+1)-fase(tau)); % salto de 10 graus na fase em t=tau
SNRdB=51; % SRN alvo em dB
x=A*cos(fase); % sinal sem ruído
vx=var(x); % variância do ruído
cn=sqrt(vx/10^(SNRdB/10))*randn(size(x)); % componente de ruído
xn=x+cn; % sinal com ruído
SNR_med=10*log10(vx/var(cn))
plot(x);
figure
plot(gradient(fase))

xa=hilbert(xn); 
phase_est=phase(xa);

figure 
fi_anomaly=gradient(phase_est);
fi_anomaly=fi_anomaly-median(fi_anomaly(10:end-10));
fi_anomaly=fi_anomaly(10:end-10);
%plot(10:length(fi_anomaly)-10,abs(fi_anomaly(10:end-10)))
plot(9+(1:length(fi_anomaly)),abs(fi_anomaly))
[v,tau_est]=max(abs(fi_anomaly));
tau_est(1)+9

