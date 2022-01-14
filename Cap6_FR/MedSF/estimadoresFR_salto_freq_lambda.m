%testes para estimação de fase, frequência e ROCOF

%variacoes interessantes:
% f1 = 59, 61, etc
% dois saltos
% encontrar um lambda otimo para o PATV

clear all; close all; clc
SNR = 60;
%fixed parameters
F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency
KaS = 0.0; %[degrees]
KxS = 0.0; % [relative magnitude step]
KfS = 1.; %[Hz] %size of frequency step
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

lambda_step = .1;
n_lambdas = 100; la_ini = .001;
lambda_n = la_ini+(0:lambda_step:(n_lambdas-1)*lambda_step);

tau_vec = 0.1:0.1:0.5;
%tau_vec = 0.5

for tau_i = 1:length(tau_vec)
for k=1:n_lambdas
for j = 1:ncurves
    Signal = SigGEN2(F0,F1,Fs,phi_n(j),NCycles,tau_vec(tau_i),tau2,SNR,KaS, KxS,KfS,nbits);
    z=hilbert(Signal');
    Psi_i = unwrap(angle(z)); %[rad]
    %Psi_i = phase(z);
    f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    az = abs(z);
      
    tau = tau_vec(tau_i)*T;
    Fref = (tau*F1 + (T-tau)*(F1 + KfS))/T;
    tau_n1 = floor(tau*NSamples/T);
    
    %estimador EFx
    %[f1_est,f2_est,F_est,fu,ri] = EF4(Psi_i,az,Fs,tau_n1,lambda_n(k));
    [f1_est,f2_est,F_est,fu,ri] = EF5(f_i,az,tau_n1,lambda_n(k));
    %[f1_est,f2_est,F_est,fu,ri] = EF6(f_i,az,tau_n1,lambda_n(k));
    
    
    ROCOF_ref = (Fref - F1)/T;
    FE(j,:) = F_est - Fref;
    FE1(j,:) = f1_est - F1;
    FE2(j,:) = f2_est - (F1 + KfS);
    kfE(j,:) = f2_est - f1_est - KfS;
    ROCOF_est = (F_est - Fref)./T;
    RFE(j,:) = ROCOF_est - ROCOF_ref; % [Hz/s]
    
end
FEmax(k,tau_i) = max(abs(FE));
FE1max(k,tau_i) = max(abs(FE1));
FE2max(k,tau_i) = max(abs(FE2));
kfEmax(k,tau_i) = max(abs(kfE));

FEmean(k,tau_i) = mean(FE);


end
end

FE
RFE
% 
% figure(1)
% semilogy(lambda_n,kfEmax,'o-')
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$KfE_{max}$ [Hz]')
% legend('\tau = 0.5T')
% legend('\tau = 0.1T','\tau = 0.2T', '\tau = 0.3T','\tau = 0.4T','\tau = 0.5T')

figure
semilogy(lambda_n,abs(FEmean),'o-')
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
xlabel('$\lambda$')
ylabel('$|\mu(FE)|$ [Hz]')
legend('\tau = 0.1T','\tau = 0.2T', '\tau = 0.3T','\tau = 0.4T','\tau = 0.5T')

figure
plot(lambda_n,abs(FEmax),'o-')
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
xlabel('$\lambda$')
ylabel('$|FE_{max}|$ [Hz]')
legend('\tau = 0.1T','\tau = 0.2T', '\tau = 0.3T','\tau = 0.4T','\tau = 0.5T')
% 
% figure
% semilogy(lambda_n,FE1max,'o-')
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$FE1_{max}$ [Hz]')
% legend('\tau = 0.1T','\tau = 0.2T', '\tau = 0.3T','\tau = 0.4T','\tau = 0.5T')
% 
% figure
% plot(lambda_n,FE2max,'o-')
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$FE2_{max}$ [Hz]')
% legend('\tau = 0.1T','\tau = 0.2T', '\tau = 0.3T','\tau = 0.4T','\tau = 0.5T')