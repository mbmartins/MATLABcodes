
%% Enhanced LoPATV 
% Lp minimization
clear all; close all; clc

fs=4800; % freq de amostragem
N=480; % tamenho do bloco
n=0:N-1; n=n(:); %base de tempo
Xm=1; % magnitude constante
f1=60;  % freq em Hz antes do salto em frequencia
hf=-1;   % altura do salto em freq
q=floor(0.05*N);

w1=2*pi*f1/fs;
w2=2*pi*hf/fs;

MC=500; % número de rodadas de Monte-Carlo
%sorteia phi_0
phi_0v=randi([0 179],MC,1);  % inteiro, por simplicidade
% sorteia taun
taunv= randi([48 432],MC,1);

% Setup PATV
Nit = 200; 
lambda = 0.5;
deg = 0;    % grau de do polinômio p 

f1_est=zeros(MC,1);
f2_est=zeros(MC,1);
hf_est=zeros(MC,1);
fr_est=zeros(MC,1);
fr_real=zeros(MC,1);
er_fr=zeros(MC,1);

for jj=1:MC, 
    phi_0_rad=phi_0v(jj)*pi/180; % phi_0 sorteado
    taun=taunv(jj);   % taun sorteado
    %         phi_0_rad = 0; taun = 48; %valores fixos p DEBUG
    PHI=w1*n+w2.*(n-taun).*(n>=taun)+phi_0_rad; % fase instantanea
    x=Xm.*cos(PHI);  % sinal x[n] 
    vx=var(x); SNR=60; vruido=vx./(10^(SNR/10)); 
    xn=x+sqrt(vruido)*randn(size(n));  % sinal AC ruidoso
    % Análise de Hilbert
    xa=hilbert(xn); ai=abs(xa); fasei=(phase(xa)); freqi=gradient(fasei); 
    freqi=fs*freqi/(2*pi);  % freq inst em Hz
   
    % comparacao
    [f1_est_M(jj),f2_est_M(jj),F_est_M(jj),fu,ri] = MedSF_PATV(freqi,ai,taun,lambda);

    freqi=freqi(q:end-q);  % ignorar amostras do inicio e final 
    % análise/decomposição PATV de freqi
    [freqi_est, p, cost] = patv_MM(freqi, deg, lambda, Nit);
       
    %freqi_estHz=fs*(freqi_est+p)./(2*pi);  % estimativa de f_i[n] em Hz
    freqi_estHz=(freqi_est+p);%./(2*pi);  % estimativa de f_i[n] em Hz

    if 0
    plot(fs*freqi./(2*pi));
    hold on
    freqi_estHz=fs*freqi_est./(2*pi); %plot(freqi_estHz,'r')
    pHz=fs*p./(2*pi);
    plot(freqi_estHz+0,'m','linewidth',3)
    plot(pHz,'k')
    pause
    close all
    end

    pHz=fs*p./(2*pi);
    
    f1_est(jj)=median(freqi_estHz(1:taun-q)); % estimativa de f1
    f2_est(jj)=median(freqi_estHz((taun-q+2):end));% estimativa de f2
    hf_est(jj)=f2_est(jj)-f1_est(jj);% estimativa de hf
    fr_est(jj)=(taun*f1_est(jj)+(N-taun)*f2_est(jj ))/N; % estimativa de fr
    fr_real(jj)=((taun*f1)+(N-taun)*(f1+hf))/N; % fr real
    er_fr(jj)=fr_est(jj)-fr_real(jj);  % erro de fr
    %er_fr(jj)=mean(pHz)-fr_real(jj);  % erro de fr
    
    er_fr_M(jj) =F_est_M(jj) - fr_real(jj); 
end

F1E= abs(mean(f1_est-f1))
F2E= abs(mean(f2_est-(f1+hf)))
hfE= abs(mean(hf_est-hf))
FE= abs(mean(er_fr))
FE_M= abs(mean(er_fr_M))

% pontos importantes que justificam as divergencias de resultados:
% 0 - encontrei um erro no meu codigo: o calculo de fref nao estava sendo
% atualizado a cada rodada de MC
% 1 - uso de Nit = 20 ou 200
% 2 - uso do PATV com valores de freq_i em Hz ou rad/s causam resposta
% diferente
% 3 - meu codigo esta retirando 24 amostras no inicio e 24 amostras no
% final, o que da um vetor de 432 amostras; enquanto no vetor freqi há 433
% amostras