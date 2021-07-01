function [Fraw,f1raw,f2raw,kfraw,riraw,draw,tau_nraw] = MC_estimation_Robustez(MCiter,F0,F1,Fs,phi_n,NCycles,tau_vec,SNR,KaS, KxS,KfS,nbits)
% phi_0 e tau1 tomados aleatoriamente de uma distribuicao uniforme a cada
% iteracao (vetores informados na entrada com todos os valores)
T = NCycles/F0;
NSamples = floor(NCycles*Fs/F0);

tau_n2 = floor(NSamples);

n = 1:NSamples;

for k = 1:MCiter
    
    tau1 = tau_vec(1,k); tau2 = 1.0; %tau_vec(2,k);
    tau_n1 = floor(tau1*NSamples);    
    tau_n2 = floor(tau2*NSamples);    
    
    Fref(k) = (tau1*T*F1 + (T - tau1*T)*(F1 + KfS))/T; %one freq step only
    ROCOF_ref = (Fref - F1)/T;

    Signal = SigGEN2(F0,F1,Fs,phi_n(k),NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits);
    z=hilbert(Signal');
    Psi_i = unwrap(angle(z)); %[rad]
    %wvd(Signal,Fs,'MinThreshold',10)  %resolução muito baixa ~5Hz
    %Psi_i = phase(z);
    f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    az = abs(z);

    % Estimation and analysis
    [f1_est1(k),f2_est1(k),F_est1(k),fu1,ri1] = EF1(f_i,tau_n1);
    [dmax1(k),tau_n_est1(k)] = max(abs(ri1));
    % valor default (ex:quando todos sao zero)
    % deixando igual ao salto de freq para commparacao
    lambdaEF5 = 0.5; %otimizado para min FE1
        
    [f1_est5(k),f2_est5(k),F_est5(k),fu5,ri5,ru5] = EF5(f_i,az,tau_n1,lambdaEF5);
    [dmax5(k),tau_n_est5(k)] = max(abs(ri5));

    kf1(k) = f2_est1(k) - f1_est1(k);
    kf5(k) = f2_est5(k) - f1_est5(k);
end

Fraw = [F_est1; F_est1; F_est1; F_est1; F_est5; F_est1] - Fref;
f1raw = [f1_est1; f1_est1; f1_est1; f1_est1; f1_est5; f1_est1] - F1;
f2raw = [f2_est1; f2_est1; f2_est1; f2_est1; f2_est5; f2_est1] - (F1 + KfS);
kfraw = [kf1; kf1; kf1; kf1; kf5; kf1] - KfS;
riraw = [ri1'; ri1'; ri1'; ri1'; ri5'; ri1'];
draw = [dmax1' dmax1' dmax1' dmax1' dmax1' dmax1'];
%th = [0.01 0.01 0.01 0.01 0.01 0.01]; 
% % criterio da deteccao do pico
% teste_d = draw > th; 
% % criterio da validacao do tamanho do salto
% salto_minimo = 0.02; %[Hz]
% teste_k = kfraw' + KfS > salto_minimo;
% tau_nraw = [tau_n_est1' tau_n_est2' tau_n_est3' tau_n_est4' tau_n_est5' tau_n_est6'].*teste_d.*teste_k;

