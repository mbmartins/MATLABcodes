%testes para estimação de fase, frequência e ROCOF

%variacoes interessantes:
% f1 = 59, 61, etc
% dois saltos
% encontrar um lambda otimo para o PATV

clear all; close all; clc
SNR = 600;
%fixed parameters
F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency
KaS = 10.0; %[degrees]
KxS = 0.0; % [relative magnitude step]
KfS = 0.0; %[Hz] %size of frequency step
Fs =4800;
%Fs = 5120*60;
NCycles = 6;
T = NCycles/F0;
phi_0 = 0; %angle phi_0 [degrees]
tau1 = 0.5;  % in proportion of T
tau2 = 1.; % in proportion of T; set tau2 = 1 if you dont want two steps
NSamples = floor(NCycles*Fs/F0);
tau_n1 = floor(tau1*NSamples);
tau_n2 = floor(tau2*NSamples);
nbits = 16;

phistep = 5; %[degrees]
ncurves = 1;
phi_n = (0:phistep:(ncurves-1)*phistep) + phi_0;

tau_vec = (0.1:0.005:0.9);
tau_n = floor(tau_vec*NSamples);

for tau_i = 1:length(tau_vec)
for j = 1:ncurves
    Signal = SigGEN2(F0,F1,Fs,phi_n(j),NCycles,tau_vec(tau_i),tau2,SNR,KaS, KxS,KfS,nbits);
    z=hilbert(Signal');
    Psi_i = unwrap(angle(z)); %[rad]
    %Psi_i = phase(z);
    f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    az = abs(z);

    % Estimation and analysis
    Fref = (tau1*T*F1 + (T - tau1*T)*(F1 + KfS))/T; %one step only
    ROCOF_ref = (Fref - F1)/T;
    n = 1:NSamples;
    br = 0.05; brn = floor(br*NSamples); % number of samples to be ignored
    %Estimator EF1 
    f_u1 = f_i;
    F_est(1) = median(f_u1);

    %Estimator EF2
    f_u2 = f_i.*az./median(az);
    f_u2 = f_u2(brn:end-brn-1);
    F_est(2) = median(f_u2);

    %Estimator EF3 - considerando dois taus somente
    f_u3 = [f_u2(1:tau_n1-2*brn); f_u2(tau_n1+1:tau_n2-2*brn); f_u2(tau_n2+1:end)];
    %f_u3 = [f_u1(1:tau_n-brn); f_u1(tau_n+1:end)];
    F_est(3) = median(f_u3);
    %F_est(3) = mean(f_u3);
    
    %Estimador EF4 - procurar um valor ótimo de lambda que minimize o erro
    %médio/máximo de Fr ??
    d = 1; lambda = 0.09; Nit = 10; %%% aplicar aqui o PATV, como no AMPS estendido
    [x, p, cost, u, v] = patv_MM(Psi_i, d, lambda, Nit);
    Psi_PATV = p;
    f_u4=gradient(Psi_PATV)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    F_est(4) = median(f_u4); 
    phi_0_est(j) = Psi_PATV(1)*180/pi;
    
    tau = tau_vec(tau_i)*T;
    Fref = (tau*F1 + (T-tau)*(F1 + KfS))/T;
    ROCOF_ref = (Fref - F1)/T;
    
    FE(tau_i,:) = F_est - Fref;
    ROCOF_est = (F_est - Fref)./T;
    RFE(tau_i,:) = ROCOF_est - ROCOF_ref; % [Hz/s]
    
end
end

FE
RFE

FE_ = mean(FE)
RFE_ = mean(RFE)

figure
plot(tau_n,FE,'.-')
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
xlabel('$\tau_n$')
ylabel('$FE$ [Hz]')
legend('EF1','EF2','EF3','EF4')
