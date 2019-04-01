clear all; close all; clc
fs=5000;% frequencia de amostragem 
f0=60; % frequencia nominal em Hz
A=10; % magnitude nominal 
fasei=90; % fase inicial em graus
phi=fasei*pi/180; % fase inicial em rad 
N=500; 
pertau=.5;  % localiza��o do salto em termos do percentual do tamanho da janela 
%tau=round(N*pertau); % �ndice do instante (arredondado) de ocorr�ncia do salto
tau=212;  % tau ajustado para o salto de fase coincidir com um maximo local do seno
n=(0:(N-1))'; % base de tempo unit�rio
fase=2*pi*f0.*n/fs+phi; % fase limpa
fase(tau:end)=fase(tau:end)+10*pi/180; % salto de 10 graus na fase em t=tau (conforme eq (2))
%fase(tau:end)=fase(tau:end)+10*pi/180-(fase(tau+1)-fase(tau)); % salto de 10 graus na fase em t=tau
SNRdB=60; % SRN alvo em dB
x=A*cos(fase); % sinal sem ru�do
if 1
vx=var(x); % vari�ncia do ru�do
cn=sqrt(vx/10^(SNRdB/10))*randn(size(x)); % componente de ru�do
xn=x+cn; % sinal com ru�do
SNR_med=10*log10(vx/var(cn))
end
%xn=x;
%plot(x);
figure
% Compensa��o do perfil de oscila��o da fi devido � presen�a do salto e
% outras n�o-idealidades da an�lise de Hilbert: usa conhecimento
% pr�vio da fase n�o-ruidosa e com o salto com localiza��o conhecida
fi_ideal=gradient(fase); 
plot(fi_ideal); 

xa=hilbert(xn); 
%phase_est=phase(xa);
phase_est=unwrap(angle(xa));
fi_anomaly_noisy=gradient(phase_est);  % fi do sinal ruidoso 


xa_nl=hilbert(x);  % caso sem ru�do
%phase_est_nl=phase(xa_nl);
phase_est_nl=unwrap(angle(xa_nl));
fi_anomaly_noiseless=gradient(phase_est_nl); % fi do sinal sem ru�do

indices_compensation=(fi_ideal>=1.5*median(fi_ideal));  % perfil do ch�o de ru�do, exceto pela regi�o do salto

fi_anomaly_noiseless(indices_compensation)=median(fi_ideal);  % template com mascara na posi��o do salto

fi_anomaly_compensated=fi_anomaly_noisy-fi_anomaly_noiseless+median(fi_ideal);

hold on
plot(fi_anomaly_compensated,'r')
plot(fi_anomaly_noisy,'k')
legend('ideal','com compensa��o','sem compensa��o')


