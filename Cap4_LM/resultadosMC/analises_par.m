% analises dos resultados

clear all; close all; clc;
load('Cap4_LM\resultadosMC\filenames_par.mat')

%analises salto de fase 1 a 6, 13 a 18
filenames(1) 
load(filenames(1)) %par1000
SF_X_std60(1,:) = max(0.01*Phasor_mag_errstdev')
SF_phi_std60(1,:) = max(Phasor_ph_errstdev')
SF_F_std60(1,:) = max(0.01*Freq_stddev')

filenames(2)
load(filenames(2)) %par0100
SF_X_std60(2,:) = max(0.01*Phasor_mag_errstdev')
SF_phi_std60(2,:) = max(Phasor_ph_errstdev')
SF_F_std60(2,:) = max(0.01*Freq_stddev')

filenames(3)
load(filenames(3)) %par0010
SF_X_std60(3,:) = max(0.01*Phasor_mag_errstdev')
SF_phi_std60(3,:) = max(Phasor_ph_errstdev')
SF_F_std60(3,:) = max(0.01*Freq_stddev')

filenames(4)
load(filenames(4)) %par0010
SF_X_std60(4,:) = max(0.01*Phasor_mag_errstdev')
SF_phi_std60(4,:) = max(Phasor_ph_errstdev')
SF_F_std60(4,:) = max(0.01*Freq_stddev')

