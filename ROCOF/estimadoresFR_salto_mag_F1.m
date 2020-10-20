%testes para estimação de fase, frequência e ROCOF

clear all; close all; clc
%fixed parameters
SNR = 60; %Signal noise ratio[dB]

%%% ajustar os valores de k para controlar os saltos
k_a = 0.0; % phase (angle) step height [degrees]
k_x = 0.1; % magnitude step height [relative]
k_f = 0.0; % frequency step height [Hz]
k_r = 0.0; % rocof step (not yet working)

F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency

Fs = 4800; % sampling frequency [Hz]
NCycles = 6; % signal number of generated nominal cycles
T = NCycles/F0; % fundamental cycles duration
NSamples = floor(NCycles*Fs/F0); % total number of signal samples
phi_0 = 90; % initial angle phi_0 [degrees]
tau1 = 0.5;  % first step location in proportion of T
tau2 = 1.; % second step location in proportion of T; set tau2 = 1 if you dont want two steps
tau_n1 = floor(tau1*NSamples); %first step sample location
tau_n2 = floor(tau2*NSamples); %2nd step sample location
nbits = 16; % number of bits for the simulated signal generator

% phistep = 10; %[degrees]
% ncurves = 36;
% phi_n = (0:phistep:(ncurves-1)*phistep) + phi_0;

tau_vec = 0.5;
%tau_vec = (0.1:0.01:0.9);
tau_n = floor(tau_vec*NSamples);

F1_vec = [58:0.1:62];

for j = 1:length(F1_vec)
    %para acompanhamento
    j
    F1 = F1_vec(j);
    %----- Loop de Monte Carlo -----
    MC_iterations = 300;
    [Fraw,f1raw,f2raw,kfraw] = MC_estimation(MC_iterations,F0,F1_vec(j),Fs,phi_0,NCycles,tau_vec,tau2,SNR,k_a, k_x,k_f,nbits);
    F_mean = mean(Fraw,2);
    %------- calculos estatisticos ----
    f1_mean(j,:) = mean(f1raw,2);
    f1_std(j,:) = std(f1raw,0,2);
    f2_mean(j,:) = mean(f2raw,2);
    f2_std(j,:) = std(f2raw,0,2);
    kf_mean(j,:) = mean(kfraw,2);
    kf_std(j,:) = std(kfraw,0,2);
    tau = tau_vec*T;
    Fref = (tau*F1 + (T-tau)*(F1 + k_f))/T; %only one step
    ROCOF_ref = (Fref - F1)/T;
    FE(j,:) = F_mean - Fref;
    FE_std(j,:) = std(Fraw,0,2);
    FE1(j,:) = f1_mean(j,:) - F1;
    KFE(j,:) = kf_mean(j,:) - k_f;
    %ROCOF_est = (F_est - F1)./T;
    %RFE(tau_i,:) = ROCOF_est - ROCOF_ref; % [Hz/s]
    

end

c = ['k','b','c','r','m','g'];
    subplot(1,3,1)
    for k = 1:6
        plot(F1_vec,FE(:,k),c(k)); hold on;
        plot(F1_vec,FE(:,k) + FE_std(:,k),[c(k),'--']);
        plot(F1_vec,FE(:,k) - FE_std(:,k),[c(k),'--']); 
    end
    title('a)')
    xlabel('x','Interpreter','latex');
    ylabel('y','Interpreter','latex');
    xlabel('$F_1 [Hz]$')
    ylabel('$FE$ [Hz]')
    %legend('EF1','EF2','EF3','EF4','EF5','EF6')
    grid on
    
    subplot(1,3,2)
    for k = 1:6
        plot(F1_vec,FE1(:,k),c(k)); hold on;
        plot(F1_vec,FE1(:,k) + f1_std(:,k),[c(k),'--']);
        plot(F1_vec,FE1(:,k) - f1_std(:,k),[c(k),'--']); 
    end
    title('b)')
    xlabel('x','Interpreter','latex');
    ylabel('y','Interpreter','latex');
    xlabel('$F_1 [Hz]$')
    ylabel('$FE1$ [Hz]')
    %legend('EF1','EF2','EF3','EF4','EF5','EF6')
    grid on
    
    subplot(1,3,3)
    for k = 1:6
        plotkP(k) = plot(F1_vec,KFE(:,k),c(k)); hold on;
        plot(F1_vec,KFE(:,k) + kf_std(:,k),[c(k),'--']);
        plot(F1_vec,KFE(:,k) - kf_std(:,k),[c(k),'--']); 
    end
    title('c)')
    xlabel('x','Interpreter','latex');
    ylabel('y','Interpreter','latex');
    xlabel('$F_1 [Hz]$')
    ylabel('$KFE1$ [Hz]')
    legend(plotkP, 'EF1','EF2','EF3','EF4','EF5','EF6')
    grid on
    
    save('salto mag\estimadoresFR_salto_mag_F1')