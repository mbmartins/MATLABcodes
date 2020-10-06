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

lambda_step = .5;
n_lambdas = 25; la_ini = 1.;
lambda_n = la_ini+(0:lambda_step:(n_lambdas-1)*lambda_step);

tau_vec = 0.05:0.1:0.5;

for tau_i = 1:length(tau_vec)
for k=1:n_lambdas
for j = 1:ncurves
    Signal = SigGEN2(F0,F1,Fs,phi_n(j),NCycles,tau_vec(tau_i),tau2,SNR,KaS, KxS,KfS,nbits);
    z=hilbert(Signal');
    Psi_i = unwrap(angle(z)); %[rad]
    %Psi_i = phase(z);
    f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    az = abs(z);

    %Estimador EF4 - procurar um valor ótimo de lambda que minimize o erro
    %médio/máximo de Fr ??
%     d = 1; lambda = lambda_n(k); Nit = 10; %%% aplicar aqui o PATV, como no AMPS estendido
%     [x, p, cost, u, v] = patv_MM(Psi_i, d, lambda, Nit);
%     Psi_PATV = p;
%     f_u4=gradient(Psi_PATV)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
%     F_est(4) = median(f_u4); 
%     phi_0_est(j) = Psi_PATV(1)*180/pi;
       
    tau = tau_vec(tau_i)*T;
    Fref = (tau*F1 + (T-tau)*(F1 + KfS))/T;
    tau_n1 = floor(tau*NSamples/T);
    
    %estimador EF5
    [f1_est5,f2_est5,F_est,fu5,ri5] = EF5(f_i,az,tau_n1,lambda_n(k));
    
    ROCOF_ref = (Fref - F1)/T;
    FE(j,:) = F_est - Fref;
    FE1(j,:) = f1_est5 - F1;
    FE2(j,:) = f2_est5 - (F1 + KfS);
    ROCOF_est = (F_est - Fref)./T;
    RFE(j,:) = ROCOF_est - ROCOF_ref; % [Hz/s]
    
end
FEmax(k,tau_i) = max(abs(FE));
FE1max(k,tau_i) = max(abs(FE1));
FE2max(k,tau_i) = max(abs(FE2));
end
end

FE
RFE

figure
plot(lambda_n,FEmax,'o-')
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
xlabel('$\lambda$')
ylabel('$FE_{max}$ [Hz]')
legend('\tau = 0.1T','\tau = 0.2T', '\tau = 0.3T','\tau = 0.4T','\tau = 0.5T')


figure
plot(lambda_n,FE1max,'o-')
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
xlabel('$\lambda$')
ylabel('$FE1_{max}$ [Hz]')
legend('\tau = 0.1T','\tau = 0.2T', '\tau = 0.3T','\tau = 0.4T','\tau = 0.5T')

figure
plot(lambda_n,FE2max,'o-')
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
xlabel('$\lambda$')
ylabel('$FE2_{max}$ [Hz]')
legend('\tau = 0.1T','\tau = 0.2T', '\tau = 0.3T','\tau = 0.4T','\tau = 0.5T')