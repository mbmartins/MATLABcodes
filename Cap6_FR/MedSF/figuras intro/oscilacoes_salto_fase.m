clear all; close all; clc;
% geracao de grafico para ilustrar as oscilacoes em fi quando há salto de
% fase
hf = 0;
NCycles = 6;
F0 = 60;
Fs = 4800;
F1 = 60;
phi_n = 0;
tau1 = 0.5;
SNR = 100;
KaS = 0; KxS = 0.1;
T = NCycles/F0;
N = floor(NCycles*Fs/F0);
tau_n1 = floor(tau1*N);

% br = 80/480; brn = floor(br*NSamples); % number of samples to be ignored
% brmask = [(brn+1:tau_n1-brn+1) (tau_n1+brn+2:NSamples-brn)];
Fref = (tau1*T*F1 + (T - tau1*T)*(F1 + hf))/T; %one step only
ROCOF_ref = (Fref - F1)/T;
n = 1:N;

k = 1;
nbits = 16;
tau2 = 10;
    x = SigGEN2(F0,F1,Fs,phi_n,NCycles,tau1,tau2,SNR,KaS, KxS,hf,nbits);
        
%         % --- gerador alternativo
%         w1 = 2*pi*F1/Fs; w2 = 2*pi*hf/Fs; Xm = 1;
%         phi_0_rad=phi_n*pi/180; % phi_0 sorteado
%         phi_0_rad = phi_0_rad + 2*pi*(F0 - F1)*tau_n1/Fs; % fator de correção para F1 fora da nominal
%         PHI=w1*n+w2.*(n-tau_n1).*(n>=tau_n1)+phi_0_rad; % fase instantanea
%         x=Xm.*cos(PHI);  % sinal x[n] 
%         %plot(x); hold on; %plot(tau_n1,Xm,'rx')
        vx=var(x);  
        vruido=vx./(10^(SNR/10)); 
        xn=x+sqrt(vruido)*randn(size(n));  % sinal AC ruidoso]
    
    
    z=hilbert(xn);
    Psi_i = unwrap(angle(z)); %[rad]
    %wvd(Signal,Fs,'MinThreshold',10)  %resolução muito baixa ~5Hz
    %Psi_i = phase(z);
    f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
n = 1:N;
s1 = subplot(211)
q = floor(0.05*N);

    p = plot(n,f_i,'k.--',n(q:N-q),f_i(q:N-q),'r.-')
    %axis([0 480 59.2 60.8])
    xlabel({'n [amostras]';'a)'});
    ylb = ylabel('$\hat{f}_i[n],\tilde{f}_i[n]$ [Hz]');
    ylb.Interpreter = 'latex';
    grid on
    lgd = legend('$\hat{f}_i[n]$','$\tilde{f}_i[n]$');
    lgd.Interpreter = 'latex';
    s2 = subplot(212)
    p2 = plot(n,x,'k.-')
    xlabel({'n [amostras]';'b)'});
    ylb = ylabel('$x[n]$ [V]');
    ylb.Interpreter = 'latex';
    grid on