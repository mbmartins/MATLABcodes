% analises dos resultados

clear all; close all; clc;
load('Cap4_LM\resultadosMC\filenames_par_tau.mat')

%analises salto de fase
filenames(1) 
load(filenames(1)) %par1000: Xm
ui = 2; % de tau_n
ui = ui/sqrt(2);

SF_X_std60(1,:) = max(0.01*Phasor_mag_errstdev')
SF_phi_std60(1,:) = max(Phasor_ph_errstdev')
SF_F_std60(1,:) = max(0.01*Freq_stddev')

    % coeficientes de sensibilidade
    c_SF_Xmag(1) = SF_X_std60(1)/(ui)
    c_SF_Xphi(1) = SF_phi_std60(1)/(ui)
    c_SF_F(1) = SF_F_std60(1)/(ui)

%analises salto mag
filenames(2) 
load(filenames(2)) %par1000: Xm
SM_X_std60(1,:) = max(0.01*Phasor_mag_errstdev');
SM_phi_std60(1,:) = max(Phasor_ph_errstdev')
SM_F_std60(1,:) = max(0.01*Freq_stddev')

    % coeficientes de sensibilidade
    c_SM_Xmag(1) = SM_X_std60(1)/(ui)
    c_SM_Xphi(1) = SM_phi_std60(1)/(ui)
    c_SM_F_Xm(1) = SM_F_std60(1)/(ui)

