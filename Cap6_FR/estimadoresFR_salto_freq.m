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
KaS = 0.0; %[degrees]
KxS = 0.0; % [relative magnitude step]
KfS = 1.0; %[Hz] %size of frequency step
Fs =4800;
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
ncurves = 180/phistep;
phi_n = (0:phistep:(ncurves-1)*phistep) + phi_0;

for j = 1:ncurves
    Signal = SigGEN2(F0,F1,Fs,phi_n(j),NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits);
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
    d = 1; lambda = 5; Nit = 10; %%% aplicar aqui o PATV, como no AMPS estendido
    [x, p, cost, u, v] = patv_MM(Psi_i, d, lambda, Nit);
    Psi_PATV = p;
    f_u4=gradient(Psi_PATV)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    F_est(4) = median(f_u4); 
    phi_0_est(j) = Psi_PATV(1)*180/pi
    
    FE(j,:) = F_est - Fref;
    ROCOF_est = (F_est - Fref)./T;
    RFE(j,:) = ROCOF_est - ROCOF_ref; % [Hz/s]
end
FE
RFE

figure
plot(phi_n,(FE-mean(FE))/F1)
xlabel('\phi_0 [graus]')
ylabel('FE - \mu(FE) [Hz/Hz]')
legend('EF1','EF2','EF3','EF4')

figure
plot(phi_n,FE)
xlabel('\phi_0 [graus]')
ylabel('FE [Hz]')
legend('EF1','EF2','EF3','EF4')
% 
% figure
% plot(phi_n,phi_0_est)
% xlabel('\phi_0 [graus]')
% ylabel('phi_0_est [graus]')

figure
fref = Fref*ones(size(f_u1));
plot(f_u1); hold on; plot(f_u2); plot(f_u3); plot(f_u4);plot(fref,'k')
legend('fu1','fu2','fu3','fu4')