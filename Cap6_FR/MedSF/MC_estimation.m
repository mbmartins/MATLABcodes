function [Fraw,f1raw,f2raw,kfraw] = MC_estimation(MC_iterations,F0,F1,Fs,phi_n,NCycles,tau1,tau2,SNR,KaS, KxS,hf,nbits)
% OBS: SOMENTE PARA SALTOS DE FREQUENCIA

T = NCycles/F0;
NSamples = floor(NCycles*Fs/F0);
tau_n1 = floor(tau1*NSamples);
tau_n2 = floor(tau2*NSamples);
% br = 80/480; brn = floor(br*NSamples); % number of samples to be ignored
% brmask = [(brn+1:tau_n1-brn+1) (tau_n1+brn+2:NSamples-brn)];
Fref = (tau1*T*F1 + (T - tau1*T)*(F1 + hf))/T; %one step only
ROCOF_ref = (Fref - F1)/T;
n = 1:NSamples;

for k = 1:MC_iterations
    %Signal = SigGEN2(F0,F1,Fs,phi_n,NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits);
        
%         % --- gerador alternativo
        w1 = 2*pi*F1/Fs; w2 = 2*pi*hf/Fs; Xm = 1;
        phi_0_rad=phi_n*pi/180; % phi_0 sorteado
        phi_0_rad = phi_0_rad + 2*pi*(F0 - F1)*tau_n1/Fs; % fator de correção para F1 fora da nominal
        PHI=w1*n+w2.*(n-tau_n1).*(n>=tau_n1)+phi_0_rad; % fase instantanea
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
    az = abs(z);

    % Estimation and analysis
    [f1_est1(k),f2_est1(k),F_est1(k),fu1,ri1] = EF1(f_i,tau_n1);
    
        % salto de frequencia
        lambda = 0.5; %otimizado para min hfE
    
    [f1_est5(k),f2_est5(k),F_est5(k),fu5,ri5,ru5] = MedSF_PATV(f_i,az,tau_n1,lambda);
    
    kf1(k) = f2_est1(k) - f1_est1(k);
    kf5(k) = f2_est5(k) - f1_est5(k);
end

Fraw = [F_est1; F_est1; F_est1; F_est1; F_est5; F_est1];
f1raw = [f1_est1; f1_est1; f1_est1; f1_est1; f1_est5; f1_est1];
f2raw = [f2_est1; f2_est1; f2_est1; f2_est1; f2_est5; f2_est1];
kfraw = [kf1; kf1; kf1; kf1; kf5; kf1];