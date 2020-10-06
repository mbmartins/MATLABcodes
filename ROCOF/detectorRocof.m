%testes para estimação de fase, frequência e ROCOF
clear all; close all; clc
SNR = 600;
%angle phi_0 [degrees]
phi_0 = 0;
%fixed parameters
F0 = 60;
F1 = 60;
tau1 = 0.5;  % in proportion of the time window
tau2 = 2; % in proportion of the time window
KaS = 0.0; %[degrees]
KxS = 0.0; % [relative magnitude step]
KfS = 0.1; %[Hz] %size of frequency step
Fs =4800;
NCycles = 6;
T = NCycles/F0;
NSamples = floor(NCycles*Fs/F0);
tau_n = floor(tau1*NSamples)

phistep = 5; %[degrees]
ncurves = 360/phistep;
phi_n = (0:phistep:(ncurves-1)*phistep) + phi_0;
n = 1:NSamples;
br = 0.05; brn = br*NSamples; %indexes of samples to ignore
br_mask = floor((br)*NSamples):floor((1-br)*NSamples);

for j = 1:ncurves
    Signal(j,:) = SigGEN2(F0,F1,Fs,phi_n(j),NCycles,tau1,tau2,SNR,KaS, KxS,KfS);
    z=hilbert(Signal(j,:)');
    Psi_i(:,j) = unwrap(angle(z)); %[rad]
    %Psi_i = phase(z);
    dfh(:,j)=gradient(Psi_i(:,j))*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    mz = abs(z);
    df(:,j) = dfh(:,j).*mz./median(mz); %compensation of magnitude interference
    fest(j) = median(dfh(:,j));
    fest_comp(j) = median(df(:,j));
    ROCOF_ih(:,j) = gradient(dfh(:,j))/T;
    ROCOF_i(:,j) = gradient(df(:,j))/T;
    ROCOF_est(j) = median(ROCOF_ih(:,j));
    ROCOF_est_comp(j) = median(ROCOF_i(:,j))
    f1est(j) = median(df(1:tau_n,j));
    f2est(j) = median(df(tau_n+1:NSamples,j));   
    fest_compDF2(j) = (f2est(j) + f1est(j))/2;
    ROCOF_estDF2(j) = (f2est(j) - f1est(j))/T;

    [ROCOF_peak,tau_nn] = max(abs(ROCOF_i(br_mask)));
    tau_n_est(j) = tau_nn + brn;
end

