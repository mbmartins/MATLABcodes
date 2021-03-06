function [Fraw,f1raw,f2raw,kfraw] = MC_estimation(MC_iterations,F0,F1,Fs,phi_n,NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits)
% OBS: SOMENTE PARA SALTOS DE FREQUENCIA

T = NCycles/F0;
NSamples = floor(NCycles*Fs/F0);
tau_n1 = floor(tau1*NSamples);
tau_n2 = floor(tau2*NSamples);
% br = 80/480; brn = floor(br*NSamples); % number of samples to be ignored
% brmask = [(brn+1:tau_n1-brn+1) (tau_n1+brn+2:NSamples-brn)];
Fref = (tau1*T*F1 + (T - tau1*T)*(F1 + KfS))/T; %one step only
ROCOF_ref = (Fref - F1)/T;
n = 1:NSamples;

for k = 1:MC_iterations
    Signal = SigGEN2(F0,F1,Fs,phi_n,NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits);
    z=hilbert(Signal');
    Psi_i = unwrap(angle(z)); %[rad]
    %wvd(Signal,Fs,'MinThreshold',10)  %resolu��o muito baixa ~5Hz
    %Psi_i = phase(z);
    f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    az = abs(z);

    % Estimation and analysis
    [f1_est1(k),f2_est1(k),F_est1(k),fu1,ri1] = EF1(f_i,tau_n1);
    %[f1_est2(k),f2_est2(k),F_est2(k),fu2,ri2] = EF2(f_i,az,tau_n1);
    %[f1_est3(k),f2_est3(k),F_est3(k),fu3,ri3] = EF3(f_i,az,tau_n1);
    
        % salto de frequencia
        lambdaEF4 = 6.5; %otimizado para minimizar FE
        lambdaEF5 = 0.5; %otimizado para min FE1
        lambdaEF6 = 1.5; %otimizado para min FE1    
    
%     if abs(KxS) > 0
%     %salto de magnitude
%     lambdaEF4 = 0.3; %otimizado para minimizar FE
%     lambdaEF5 = 0.5; %otimizado para min FE1
%     lambdaEF6 = 0.5; %otimizado para min FE1    
%     end
% 
%     if abs(KaS) > 0
%     %salto de fase
%     lambdaEF4 = 0.05; %otimizado para minimizar FE
%     lambdaEF5 = 2; %otimizado para min FE1
%     lambdaEF6 = 2.0; %otimizado para min FE1    
%     end
    
    
    %[f1_est4(k),f2_est4(k),F_est4(k),fu4,ri4,ru4] = EF4(Psi_i,az,Fs,tau_n1,lambdaEF4);
    [f1_est5(k),f2_est5(k),F_est5(k),fu5,ri5,ru5] = EF5(f_i,az,tau_n1,lambdaEF5);
    %[f1_est6(k),f2_est6(k),F_est6(k),fu6,ri6,ru6] = EF6(f_i,az,tau_n1,lambdaEF6);
    
    kf1(k) = f2_est1(k) - f1_est1(k);
    %kf2(k) = f2_est2(k) - f1_est2(k);
    %kf3(k) = f2_est3(k) - f1_est3(k);
    %kf4(k) = f2_est4(k) - f1_est4(k);
    kf5(k) = f2_est5(k) - f1_est5(k);
    %kf6(k) = f2_est6(k) - f1_est6(k);
end

Fraw = [F_est1; F_est1; F_est1; F_est1; F_est5; F_est1];
f1raw = [f1_est1; f1_est1; f1_est1; f1_est1; f1_est5; f1_est1];
f2raw = [f2_est1; f2_est1; f2_est1; f2_est1; f2_est5; f2_est1];
kfraw = [kf1; kf1; kf1; kf1; kf5; kf1];