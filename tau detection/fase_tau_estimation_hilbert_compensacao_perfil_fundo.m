clear all; close all; clc
fs=5000;% frequencia de amostragem 
f0=60; % frequencia nominal em Hz
A=10; % magnitude nominal 
fasei=90; % fase inicial em graus
phi=fasei*pi/180; % fase inicial em rad 
N=500; 
pertau=.5;  % localização do salto em termos do percentual do tamanho da janela 
%tau=round(N*pertau); % índice do instante (arredondado) de ocorrência do salto
tau=212;  % tau ajustado para o salto de fase coincidir com um maximo local do seno
n=(0:(N-1))'; % base de tempo unitário
fase=2*pi*f0.*n/fs+phi; % fase limpa
fase(tau:end)=fase(tau:end)+10*pi/180; % salto de 10 graus na fase em t=tau (conforme eq (2))
%fase(tau:end)=fase(tau:end)+10*pi/180-(fase(tau+1)-fase(tau)); % salto de 10 graus na fase em t=tau
SNRdB=60; % SRN alvo em dB
x=A*cos(fase); % sinal sem ruído
if 1
vx=var(x); % variância do ruído
cn=sqrt(vx/10^(SNRdB/10))*randn(size(x)); % componente de ruído
xn=x+cn; % sinal com ruído
SNR_med=10*log10(vx/var(cn))
end
%xn=x;
%plot(x);
figure
% Compensação do perfil de oscilação da fi devido à presença do salto e
% outras não-idealidades da análise de Hilbert: usa conhecimento
% prévio da fase não-ruidosa e com o salto com localização conhecida
fi_ideal=gradient(fase); 
plot(fi_ideal); 

xa=hilbert(xn); 
%phase_est=phase(xa);
phase_est=unwrap(angle(xa));
fi_anomaly_noisy=gradient(phase_est);  % fi do sinal ruidoso 


xa_nl=hilbert(x);  % caso sem ruído
%phase_est_nl=phase(xa_nl);
phase_est_nl=unwrap(angle(xa_nl));
fi_anomaly_noiseless=gradient(phase_est_nl); % fi do sinal sem ruído

indices_compensation=(fi_ideal>=1.5*median(fi_ideal));  % perfil do chão de ruído, exceto pela região do salto

fi_anomaly_noiseless(indices_compensation)=median(fi_ideal);  % template com mascara na posição do salto

fi_anomaly_compensated=fi_anomaly_noisy-fi_anomaly_noiseless+median(fi_ideal);

hold on
plot(fi_anomaly_compensated,'r')
plot(fi_anomaly_noisy,'k')
legend('ideal','com compensação','sem compensação')


