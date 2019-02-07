clear all; close all; clc
fs=5000;% frequencia de amostragem 
f0=60; % frequencia nominal em Hz
A=10; % magnitude nominal 
fasei=90; % fase inicial em graus
phi=fasei*pi/180; % fase inicial em rad 
N=500; 
pertau=.5;  % localização do salto em termos do percentual do tamanho da janela 
%tau=round(N*pertau); % índice do instante (arredondado) de ocorrência do salto
tau=292;  % tau ajustado para o salto de fase coincidir com um maximo local do seno
n=(0:(N-1))'; % base de tempo unitário
fase=2*pi*f0.*n/fs+phi; % fase limpa
%fase(tau:end)=fase(tau:end)+10*pi/180; % salto de 10 graus na fase em t=tau (conforme eq (2))
%fase(tau:end)=fase(tau:end)+10*pi/180-(fase(tau+1)-fase(tau)); % salto de 10 graus na fase em t=tau
SNRdB=51; % SRN alvo em dB
mag=[A*ones(tau-1,1); 1.1*A*ones(N-tau+1,1)]; 
x=mag.*cos(fase); % sinal sem ruído
vx=var(x); % variância do ruído
cn=sqrt(vx/10^(SNRdB/10))*randn(size(x)); % componente de ruído
xn=x+cn; % sinal com ruído
SNR_med=10*log10(vx/var(cn))

plot(n,x,tau,x(tau),'x');
title('Signal')
figure
plot(n,gradient(fase))
xlabel('Samples');ylabel('gradient fase')

xa=hilbert(xn); 

% detector 1
%phase_est=phase(xa);
phase_est=unwrap(angle(xa));

 
fi_anomaly=gradient(phase_est);
fi_anomaly=fi_anomaly-median(fi_anomaly(10:end-10));
fi_anomaly=fi_anomaly(10:end-10);
limiar1=7*median(abs(fi_anomaly));

%plot(10:length(fi_anomaly)-10,abs(fi_anomaly(10:end-10)))
figure
plot(9+(1:length(fi_anomaly)),abs(fi_anomaly),9+(1:length(fi_anomaly)),limiar1*ones(length(fi_anomaly),1),'k')
[v,tau_est]=max(abs(fi_anomaly));
if v>limiar1,
    tau_est(1)+9
else
    disp('estimativa inválida');
end

% detector2
mag_est=abs(xa);
mag_anomaly=gradient(mag_est);
mag_anomaly=mag_anomaly-median(mag_anomaly(10:end-10));
mag_anomaly=mag_anomaly(10:end-10);
limiar2=7*median(abs(mag_anomaly));

%plot(10:length(fi_anomaly)-10,abs(fi_anomaly(10:end-10)))
figure
plot(9+(1:length(mag_anomaly)),abs(mag_anomaly),9+(1:length(mag_anomaly)),limiar2*ones(length(mag_anomaly),1),'k')
[v,tau_est]=max(abs(mag_anomaly));
if v>limiar2,
    tau_est(1)+9
else
    disp('estimativa inválida');
end

