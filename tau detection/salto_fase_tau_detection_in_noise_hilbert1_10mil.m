clear all; close all; clc

%verificar influencia:
% fase inicial
% fs
% SNR
% f0?
% pertau
% tamanho do salto

fs=5000;% frequencia de amostragem 
f0=60.0; % frequencia nominal em Hz
A=10; % magnitude nominal 
fasei=90; % fase inicial em graus
phi=fasei*pi/180; % fase inicial em rad 
ka = 10; %salto de fase em graus
N=500; 
pertau=.5;  % localiza��o do salto em termos do percentual do tamanho da janela 
%tau=round(N*pertau); % �ndice do instante (arredondado) de ocorr�ncia do salto
tau=311;  % tau ajustado para o salto de fase coincidir com um maximo local do seno
%tau=272; % tau ajustado para o salto de fase coincidir com um minimo local do seno
n=(0:(N-1))'; % base de tempo unit�rio
fase=2*pi*f0.*n/fs+phi; % fase limpa
fase(tau:end)=fase(tau:end)+ka*pi/180; % salto de ka graus na fase em t=tau (conforme eq (2))
%fase(tau:end)=fase(tau:end)+10*pi/180-(fase(tau+1)-fase(tau)); % salto de 10 graus na fase em t=tau
SNRdB=51; % SRN alvo em dB
x=A*cos(fase); % sinal sem ru�do
vx=var(x); % vari�ncia do ru�do
cn=sqrt(vx/10^(SNRdB/10))*randn(length(x),10000); % componentes de ru�do
for k=1:10000,
xn=x+cn(:,k); % sinal com ru�do
xa=hilbert(xn); 
%phase_est=phase(xa);
phase_est=unwrap(angle(xa));
fi_anomaly=gradient(phase_est);
fi_anomaly=fi_anomaly-median(fi_anomaly(10:end-10));
fi_anomaly=fi_anomaly(10:end-10);
[v,tau_est]=max(abs(fi_anomaly));
tau_est_v(k)=tau_est(1)+9;
end

%hist (tau_est_v,tau-10:tau+10)
tau_error = tau - tau_est_v;
hist(tau_error,tau_error-5:tau_error+5)
