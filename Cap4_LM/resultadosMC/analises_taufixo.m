% analises dos resultados

clear all; close all; clc;
load('Cap4_LM\resultadosMC\filenames.mat')

%analises salto de fase 1 a 6, 
%analises salto de mag 7 a 12, 

% magnitude
for j = 1:6
    filenames(j)
    load(filenames(j)) %sf 1 a 6
%     SF_X_medio93(j,:) = 0.01*Phasor_mag_errmean';
%     SF_X_max93(j,:) = 0.01*Phasor_mag_errstdev';
%     SF_phi_medio93(j,:) = Phasor_ph_errmean';
%     SF_phi_max93(j,:) = Phasor_ph_errstdev';
    SF_F_max93(j,:) = 0.01*Freq_errmax_abs';
    
    filenames(j+6)
    load(filenames(j+6)) %sf 13 a 18
%     SM_X_medio93(j,:) = 0.01*Phasor_mag_errmean';
%     SM_X_max93(j,:) = 0.01*Phasor_mag_errstdev';
%     SM_phi_medio93(j,:) = Phasor_ph_errmean';
%     SM_phi_max93(j,:) = Phasor_ph_errstdev';
%     SM_F_medio93(j,:) = 0.01*Freq_errmed';
    SM_F_max93(j,:) = 0.01*Freq_errmax_abs';     
end

keepvars = {'SF_F_medio60', 'SF_F_medio93','SF_F_std60','SF_F_max93',...
    'SF_X_medio60', 'SF_X_medio93','SF_X_std60','SF_X_max93',...
'SF_phi_medio60', 'SF_phi_medio93','SF_phi_std60','SF_phi_max93',...
'SM_F_medio60', 'SM_F_medio93','SM_F_std60','SM_F_max93',...
    'SM_X_medio60', 'SM_X_medio93','SM_X_std60','SM_X_max93',...
'SM_phi_medio60', 'SM_phi_medio93','SM_phi_std60','SM_phi_max93',...
'filenames'};
clearvars('-except',keepvars{:})

F1 = 50;
ti = 5;
tau = ti/10;

SM_F_max = max(SM_F_max93(:,ti))*F1 %[Hz]
SF_F_max = max(SF_F_max93(:,ti))*F1 %[Hz]
