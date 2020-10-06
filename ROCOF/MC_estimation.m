function [Fraw,f1raw,f2raw,kfraw] = MC_estimation(MC_iterations,F0,F1,Fs,phi_n,NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits)
T = NCycles/F0;
NSamples = floor(NCycles*Fs/F0);
tau_n1 = floor(tau1*NSamples);
tau_n2 = floor(tau2*NSamples);
br = 80/480; brn = floor(br*NSamples); % number of samples to be ignored
brmask = [(brn+1:tau_n1-brn+1) (tau_n1+brn+2:NSamples-brn)];
Fref = (tau1*T*F1 + (T - tau1*T)*(F1 + KfS))/T; %one step only
ROCOF_ref = (Fref - F1)/T;
n = 1:NSamples;

for k = 1:MC_iterations
    Signal = SigGEN2(F0,F1,Fs,phi_n,NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits);
    z=hilbert(Signal');
    Psi_i = unwrap(angle(z)); %[rad]
    %wvd(Signal,Fs,'MinThreshold',10)  %resolução muito baixa ~5Hz
    %Psi_i = phase(z);
    f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    az = abs(z);

    % Estimation and analysis
    [f1_est1(k),f2_est1(k),F_est1(k),fu1,ri1] = EF1(f_i,tau_n1);
    [f1_est2(k),f2_est2(k),F_est2(k),fu2,ri2] = EF2(f_i,az,tau_n1);
    
    %lambda = 0.36;  %otimizado para kf
    [f1_est3(k),f2_est3(k),F_est3(k),fu3,ri3] = EF3(f_i,az,tau_n1);
    
    [f1_est4(k),f2_est4(k),F_est4(k),fu4,ri4] = EF4(Psi_i,az,Fs,tau_n1);
    lambda = 6.5; %para comparar com o EF4 na mesma base
    [f1_est5(k),f2_est5(k),F_est5(k),fu5,ri5] = EF5(f_i,az,tau_n1,lambda);
    %[f1_est6(k),f2_est6(k),F_est6(k),fu6,ri6] = EF6(f_i,az,tau_n1);
    
    kf1(k) = f2_est1(k) - f1_est1(k);
    kf2(k) = f2_est2(k) - f1_est2(k);
    kf3(k) = f2_est3(k) - f1_est3(k);
    kf4(k) = f2_est4(k) - f1_est4(k);
    kf5(k) = f2_est5(k) - f1_est5(k);
%    kf6(k) = f2_est6(k) - f1_est6(k);
end

Fraw = [F_est1; F_est2; F_est3; F_est4; F_est5];
f1raw = [f1_est1; f1_est2; f1_est3; f1_est4; f1_est5];
f2raw = [f2_est1; f2_est2; f2_est3; f2_est4; f2_est5];
kfraw = [kf1; kf2; kf3; kf4; kf5];