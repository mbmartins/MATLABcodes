clear all; close all; clc;

hf = -1;
NCycles = 6;
F0 = 60;
Fs = 4800;
F1 = 60;
phi_n = 0;
phi_n2 = 45;
phi_n3 = 75;

tau1 = 0.5;
SNR = 100;
KaS = 0; KxS = 0;
T = NCycles/F0;
NSamples = floor(NCycles*Fs/F0);
tau_n1 = floor(tau1*NSamples);

% br = 80/480; brn = floor(br*NSamples); % number of samples to be ignored
% brmask = [(brn+1:tau_n1-brn+1) (tau_n1+brn+2:NSamples-brn)];
Fref = (tau1*T*F1 + (T - tau1*T)*(F1 + hf))/T; %one step only
ROCOF_ref = (Fref - F1)/T;
n = 1:NSamples;
q = floor(0.05*NSamples)
  
%         % --- gerador alternativo
        w1 = 2*pi*F1/Fs; w2 = 2*pi*hf/Fs; Xm = 1;
        phi_0_rad=phi_n*pi/180; % phi_0 sorteado
        phi_0_rad = phi_0_rad + 2*pi*(F0 - F1)*tau_n1/Fs; % fator de correção para F1 fora da nominal
        phi_0_rad2=phi_n2*pi/180; % phi_0 sorteado
        phi_0_rad2 = phi_0_rad2 + 2*pi*(F0 - F1)*tau_n1/Fs; % fator de correção para F1 fora da nominal
        phi_0_rad3=phi_n3*pi/180; % phi_0 sorteado
        phi_0_rad3 = phi_0_rad3 + 2*pi*(F0 - F1)*tau_n1/Fs; % fator de correção para F1 fora da nominal
        
        PHI=w1*n+w2.*(n-tau_n1).*(n>=tau_n1)+phi_0_rad; % fase instantanea
        x=Xm.*cos(PHI);  % sinal x[n] 
        %plot(x); hold on; %plot(tau_n1,Xm,'rx')
        vx=var(x);  
        vruido=vx./(10^(SNR/10)); 
        xn=x+sqrt(vruido)*randn(size(n));  % sinal AC ruidoso]
        
        PHI2=w1*n+w2.*(n-tau_n1).*(n>=tau_n1)+phi_0_rad2; % fase instantanea
        x2=Xm.*cos(PHI2);  % sinal x[n] 
        %plot(x); hold on; %plot(tau_n1,Xm,'rx')
        vx=var(x2);  
        vruido=vx./(10^(SNR/10)); 
        xn2=x2+sqrt(vruido)*randn(size(n));  % sinal AC ruidoso]
        
        PHI3=w1*n+w2.*(n-tau_n1).*(n>=tau_n1)+phi_0_rad3; % fase instantanea
        x3=Xm.*cos(PHI3);  % sinal x[n] 
        %plot(x); hold on; %plot(tau_n1,Xm,'rx')
        vx=var(x3);  
        vruido=vx./(10^(SNR/10)); 
        xn3=x3+sqrt(vruido)*randn(size(n));  % sinal AC ruidoso]
    
    
    z=hilbert(xn);
    z2=hilbert(xn2);
    z3=hilbert(xn3);
    
    Psi_i = unwrap(angle(z)); %[rad]
        Psi_i2 = unwrap(angle(z2)); %[rad]
            Psi_i3 = unwrap(angle(z3)); %[rad]
    
    %wvd(Signal,Fs,'MinThreshold',10)  %resolução muito baixa ~5Hz
    %Psi_i = phase(z);
    
    f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    f_i2=gradient(Psi_i2)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    f_i3=gradient(Psi_i3)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    
    n = 1:NSamples;

    ntrunc = q:NSamples - q;
n1 = q:tau_n1-1;
n2 = tau_n1+1:NSamples-q;
    f1est = median(f_i(n1));
    f2est = median(f_i(n2));
    f1estb = median(f_i3(n1));
    f2estb = median(f_i3(n2));
    
 s1 = subplot(121)
    p = plot(ntrunc,f_i(ntrunc),'b'); hold on;
    %plot(ntrunc,f_i2(ntrunc)); hold on;
    plot(ntrunc,f_i3(ntrunc),'r'); hold on;
    %axis([0 480 59.2 60.8])
    xlabel({'Amostras';'a)'});
    ylb = ylabel('$f_i[n]$');
    ylb.Interpreter = 'latex';
    grid on
    pf1a = line([q,tau_n1-1],[f1est,f1est]); pf1a.Color = 'blue'; pf1a.LineStyle = '--';
    pf1b = line([q,tau_n1-1],[f1estb,f1estb]); pf1b.Color = 'red'; pf1b.LineStyle = '--';
    
    pf2a = line([tau_n1+1,NSamples-q],[f2est,f2est]); pf2a.Color = 'blue'; pf2a.LineStyle = '--';
    pf2b = line([tau_n1+1,NSamples-q],[f2estb,f2estb]); pf2b.Color = 'red'; pf2b.LineStyle = '--';

    lg = legend('$\hat{f}_i[n], \phi_0 = 0^o$','$\hat{f}_i[n], \phi_0 = 75^o$', 'MedSF $\hat{f}_1,\hat{f}_2, \phi_0 = 0^o$', 'MedSF $\hat{f}_1, \hat{f}_2, \phi_0 = 75^o$')
    lg.Interpreter = 'latex';
    
    f2est = median(f_i(n2));
    
    s2 = subplot(122)
    h = histogram(f_i(ntrunc)); hold on;
    h2 = histogram(f_i3(ntrunc))
    h.Orientation = 'horizontal'; h2.Orientation = 'horizontal';
    h.BinLimits = s1.YLim; h2.BinLimits = s1.YLim;
    s2.YLim = h.BinLimits; s2.YLim = h2.BinLimits;
    h.NumBins = 90;  h2.NumBins = h.NumBins; 
    xlabel({'Ocorrências';'b)'});
    ylb = ylabel('$f_i[n]$');
    ylb.Interpreter = 'latex';
    grid on