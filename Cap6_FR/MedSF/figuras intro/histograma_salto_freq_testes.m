clear all; close all; clc;

hf = -1;
NCycles = 6;
F0 = 60;
Fs = 4800;
F1 = 60;
phi_n = 0;
tau1 = 0.5;
SNR = 100;
KaS = 0; KxS = 0;
T = NCycles/F0;
N = floor(NCycles*Fs/F0);
tau_n = floor(tau1*N);

% br = 80/480; brn = floor(br*NSamples); % number of samples to be ignored
% brmask = [(brn+1:tau_n1-brn+1) (tau_n1+brn+2:NSamples-brn)];
Fref = (tau1*T*F1 + (T - tau1*T)*(F1 + hf))/T; %one step only
ROCOF_ref = (Fref - F1)/T;
n = 1:N;

k = 1
    %Signal = SigGEN2(F0,F1,Fs,phi_n,NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits);
        
%         % --- gerador alternativo
        w1 = 2*pi*F1/Fs; w2 = 2*pi*hf/Fs; Xm = 1;
        phi_0_rad=phi_n*pi/180; % phi_0 sorteado
        phi_0_rad = phi_0_rad + 2*pi*(F0 - F1)*tau_n/Fs; % fator de correção para F1 fora da nominal
        PHI=w1*n+w2.*(n-tau_n).*(n>=tau_n)+phi_0_rad; % fase instantanea
        x=Xm.*cos(PHI);  % sinal x[n] 
        %plot(x); hold on; %plot(tau_n1,Xm,'rx')
        vx=var(x);  
        vruido=vx./(10^(SNR/10)); 
        xn=x+sqrt(vruido)*randn(size(n));  % sinal AC ruidoso]
    
    
    z=hilbert(xn);
    Psi_i = unwrap(angle(z)); %[rad]
    %wvd(Signal,Fs,'MinThreshold',10)  %resolução muito baixa ~5Hz
    %Psi_i = phase(z);
    f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
n = 1:N;
q = floor(0.05*N);
eps_tau = 10;

s1 = subplot(121)
    pa = plot(n(q:tau_n + eps_tau),f_i(q:tau_n+eps_tau) - F1,'b.-'); hold on;
    pa2 = plot(n(q:tau_n - eps_tau),f_i(q:tau_n-eps_tau) - F1,'r.-'); hold on;
    plot(tau_n,f_i(tau_n) - F1,'kx');
    lgd = legend('$\epsilon_\tau = -10$','$\epsilon_\tau = +10$','$\tilde{f}_i[\tau_n]$')
    lgd.Interpreter= 'latex';
    %axis([0 480 59.2 60.8])
    xlabel({'n';'a)'});
    ylb = ylabel('$\tilde{f}_i[n] - f_1$');
    ylb.Interpreter = 'latex';
    grid on
    s2 = subplot(122)
    pb = plot(n(tau_n - eps_tau:N-q),f_i(tau_n-eps_tau:N-q) - F1 - hf,'b.-'); hold on
    pb2 = plot(n(tau_n + eps_tau:N-q),f_i(tau_n + eps_tau:N-q) - F1 - hf,'r.-'); hold on;
    plot(tau_n,f_i(tau_n) - F1 - hf,'kx');
    lgd2 = legend('$\epsilon_\tau = -10$','$\epsilon_\tau = +10$','$\tilde{f}_i[\tau_n]$')
    lgd2.Interpreter= 'latex';    
%     h = histogram(f_i)
%     h.Orientation = 'horizontal';
%     h.BinLimits = s1.YLim;
%     s2.YLim = h.BinLimits;
%     h.NumBins = 45;
%     xlabel({'Ocorrências';'b)'});
    ylb = ylabel('$\tilde{f}_i[n] - f_2$');
    ylb.Interpreter = 'latex';
    grid on
        xlabel({'n';'b)'});