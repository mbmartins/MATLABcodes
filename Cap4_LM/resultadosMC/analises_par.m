% analises dos resultados

clear all; close all; clc;
load('Cap4_LM\resultadosMC\filenames_par.mat')

ui = 1/100;

%analises salto de fase
filenames(1) 
load(filenames(1)) %par1000: Xm
SF_X_std60(1,:) = max(0.01*Phasor_mag_errstdev')
SF_phi_std60(1,:) = max(Phasor_ph_errstdev')
SF_F_std60(1,:) = max(0.01*Freq_stddev')

    % coeficientes de sensibilidade
    c_SF_Xmag(1) = SF_X_std60(1)/(ui)
    c_SF_Xphi(1) = SF_phi_std60(1)/(ui)
    c_SF_F(1) = SF_F_std60(1)/(ui)

filenames(2)
load(filenames(2)) %par0100: w
SF_X_std60(2,:) = max(0.01*Phasor_mag_errstdev')
SF_phi_std60(2,:) = max(Phasor_ph_errstdev')
SF_F_std60(2,:) = max(0.01*Freq_stddev')

    c_SF_Xmag(2)= SF_X_std60(2)/(ui)
    c_SF_Xphi(2)= SF_phi_std60(2)/(ui)
    c_SF_F(2)= SF_F_std60(2)/(ui)

filenames(3)
load(filenames(3)) %par0010: phi0
SF_X_std60(3,:) = max(0.01*Phasor_mag_errstdev')
SF_phi_std60(3,:) = max(Phasor_ph_errstdev')
SF_F_std60(3,:) = max(0.01*Freq_stddev')

    c_SF_Xmag(3) = SF_X_std60(3)/(ui)
    c_SF_Xphi(3) = SF_phi_std60(3)/(ui)
    c_SF_F(3) = SF_F_std60(3)/(ui)

filenames(4)
load(filenames(4)) %par0001: ha
SF_X_std60(4,:) = max(0.01*Phasor_mag_errstdev')
SF_phi_std60(4,:) = max(Phasor_ph_errstdev')
SF_F_std60(4,:) = max(0.01*Freq_stddev')

    c_SF_Xmag(4) = SF_X_std60(4)/(ui)
    c_SF_Xphi(4) = SF_phi_std60(4)/(ui)
    c_SF_F(4) = SF_F_std60(4)/(ui)

    
%analises salto mag
filenames(5) 
load(filenames(5)) %par1000: Xm
SM_X_std60(1,:) = max(0.01*Phasor_mag_errstdev');
SM_phi_std60(1,:) = max(Phasor_ph_errstdev')
SM_F_std60(1,:) = max(0.01*Freq_stddev')

    % coeficientes de sensibilidade
    c_SM_Xmag(1) = SM_X_std60(1)/(ui)
    c_SM_Xphi(1) = SM_phi_std60(1)/(ui)
    c_SM_F_Xm(1) = SM_F_std60(1)/(ui)

filenames(6)
load(filenames(6)) %par0100: hm
SM_X_std60(2,:) = max(0.01*Phasor_mag_errstdev')
SM_phi_std60(2,:) = max(Phasor_ph_errstdev')
SM_F_std60(2,:) = max(0.01*Freq_stddev')

    c_SM_Xmag(2)= SM_X_std60(2)/(ui)
    c_SM_Xphi(2)= SM_phi_std60(2)/(ui)
    c_SM_F(2)= SM_F_std60(2)/(ui)

filenames(7)
load(filenames(7)) %par0010: w
SM_X_std60(3,:) = max(0.01*Phasor_mag_errstdev')
SM_phi_std60(3,:) = max(Phasor_ph_errstdev')
SM_F_std60(3,:) = max(0.01*Freq_stddev')

    c_SM_Xmag(3) = SM_X_std60(3)/(ui)
    c_SM_Xphi(3) = SM_phi_std60(3)/(ui)
    c_SM_F(3) = SM_F_std60(3)/(ui)

filenames(8)
load(filenames(8)) %par0001: phi0
SM_X_std60(4,:) = max(0.01*Phasor_mag_errstdev')
SM_phi_std60(4,:) = max(Phasor_ph_errstdev')
SM_F_std60(4,:) = max(0.01*Freq_stddev')

    c_SM_Xmag(4) = SM_X_std60(4)/(ui)
    c_SM_Xphi(4) = SM_phi_std60(4)/(ui)
    c_SM_F(4) = SM_F_std60(4)/(ui)
