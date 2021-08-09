function [filename] = teste_MedSF_phi0(MC_iterations,k_f,pasta)

%testes para estimação de fase, frequência e ROCOF
% em funcao de phi0, com tau fixo
%fixed parameters
SNR = 60; %Signal noise ratio[dB]

%%% ajustar os valores de k para controlar os saltos
k_a = 0.0; % phase (angle) step height [degrees]
k_x = 0.0; % magnitude step height [relative]
%k_f = 1.0; % frequency step height [Hz]
k_r = 0.0; % rocof step (not yet working)

F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency
F2 = F1 + k_f;

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
ncurves = 180/phistep;
phi_n = (0:phistep:(ncurves-1)*phistep) + phi_0;
Fref = (tau1*T*F1 + (T - tau1*T)*(F1 + k_f))/T; %one step only
ROCOF_ref = (Fref - F1)/T;

for j = 1:ncurves
    phi_n(j)
    %MC_iterations = 300;
    [Fraw,f1raw,f2raw,kfraw] = MC_estimation(MC_iterations,F0,F1,Fs,phi_n(j),NCycles,tau1,tau2,SNR,k_a, k_x,k_f,nbits);
    F_mean = mean(Fraw,2);
    
    f1_mean(j,:) = mean(f1raw,2);
    f1_std(j,:) = std(f1raw,0,2);
    
    f2_mean(j,:) = mean(f2raw,2);
    f2_std(j,:) = std(f2raw,0,2);
    kf_mean(j,:) = mean(kfraw,2);
    kf_std(j,:) = std(kfraw,0,2);
    if MC_iterations>1
        FE_std(j,:) = std(Fraw,0,2);
    else
        FE_std(j,:) = 0;
    end
    
    FE_ruido(j,:) = F_mean - Fref;
    FE1_ruido(j,:) = f1_mean(j,:) - F1;
    FE2_ruido(j,:) = f2_mean(j,:) - F2;
    KE_ruido(j,:) = kf_mean(j,:) - k_f;
    
end

% filename = pasta + "salto_freq_phi0.mat"
% save(filename)