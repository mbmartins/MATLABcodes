%testes para estimação de fase, frequência e ROCOF

clear all; close all; clc
%fixed parameters
SNR = 60; %Signal noise ratio[dB]

%%% ajustar os valores de k para controlar os saltos
k_a = 0.0; % phase (angle) step height [degrees]
k_x = 0.0; % magnitude step height [relative]
k_f = 1.0; % frequency step height [Hz]
k_r = 0.0; % rocof step (not yet working)

%parametro para o PATV
% lambda_fase = 0.09; %otimizado para salto de fase
% lambda_mag = 0.065; %otimizado para salto de fase
% lambda_freq = 6.5; %otimizado para salto de frequencia
% lambda = lambda_fase;

F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency

Fs = 4800; % sampling frequency [Hz]
NCycles = 6; % signal number of generated nominal cycles
T = NCycles/F0; % fundamental cycles duration
NSamples = floor(NCycles*Fs/F0); % total number of signal samples
phi_0 = 0; % initial angle phi_0 [degrees]
tau1 = 0.5;  % first step location in proportion of T
tau2 = 1.; % second step location in proportion of T; set tau2 = 1 if you dont want two steps
tau_n1 = floor(tau1*NSamples); %first step sample location
tau_n2 = floor(tau2*NSamples); %2nd step sample location
nbits = 16; % number of bits for the simulated signal generator

phistep = 5; %[degrees]
ncurves = 1;
phi_n = (0:phistep:(ncurves-1)*phistep) + phi_0;

tau_vec = (0.1:0.005:0.9);
tau_n = floor(tau_vec*NSamples);

for tau_i = 1:length(tau_vec)
for j = 1:ncurves
    %----- geracao do sinal real e analitico
    Signal = SigGEN2(F0,F1,Fs,phi_n(j),NCycles,tau_vec(tau_i),tau2,SNR,k_a, k_x,k_f,nbits);
    z=hilbert(Signal');
    Psi_i = unwrap(angle(z)); %[rad]
    f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    az = abs(z);

    %----- Loop de Monte Carlo -----
    MC_iterations = 300;
    [Fraw,f1raw,f2raw,kfraw] = MC_estimation(MC_iterations,F0,F1,Fs,phi_n(j),NCycles,tau_vec(tau_i),tau2,SNR,k_a, k_x,k_f,nbits);
    F_mean = mean(Fraw,2);
    %------- calculos estatisticos ----
    f1_mean(j,:) = mean(f1raw,2);
    f1_std(j,:) = std(f1raw,0,2);
    f2_mean(j,:) = mean(f2raw,2);
    f2_std(j,:) = std(f2raw,0,2);
    kf_mean(j,:) = mean(kfraw,2);
    kf_std(j,:) = std(f2raw,0,2);
    tau = tau_vec(tau_i)*T;
    Fref = (tau*F1 + (T-tau)*(F1 + k_f))/T; %only one step
    ROCOF_ref = (Fref - F1)/T;
    FE(tau_i,:) = F_mean - Fref;
    %ROCOF_est = (F_est - F1)./T;
    %RFE(tau_i,:) = ROCOF_est - ROCOF_ref; % [Hz/s]
    
end
end

FE_ = mean(FE)
%RFE_ = mean(RFE)

figure
plot(tau_n,FE)
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
xlabel('$\tau_n$')
ylabel('$FE$ [Hz]')
legend('EF1','EF2','EF3','EF4','EF5','EF6')
