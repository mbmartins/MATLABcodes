% analises dos resultados

clear all; close all; clc;
load('Cap4_LM\resultadosMC\filenames.mat')

%analises salto de fase 1 a 6, 13 a 18



%analises salto de mag 7 a 12, 19 a 24

% magnitude
for j = 1:6
    filenames(j)
    load(filenames(j)) %sf 1 a 6
    SF_X_medio60(j,:) = 0.01*Phasor_mag_errmean';
    SF_X_std60(j,:) = 0.01*Phasor_mag_errstdev';
    SF_phi_medio60(j,:) = Phasor_ph_errmean';
    SF_phi_std60(j,:) = Phasor_ph_errstdev';
    SF_F_medio60(j,:) = 0.01*Freq_errmed';
    SF_F_std60(j,:) = 0.01*Freq_stddev';    
    
    filenames(j+6)
    load(filenames(j+6)) % sm 7 a 12
    SM_X_medio60(j,:) = 0.01*Phasor_mag_errmean';
    SM_X_std60(j,:) = 0.01*Phasor_mag_errstdev';
    SM_phi_medio60(j,:) = Phasor_ph_errmean';
    SM_phi_std60(j,:) = Phasor_ph_errstdev';
    SM_F_medio60(j,:) = 0.01*Freq_errmed';
    SM_F_std60(j,:) = 0.01*Freq_stddev';

    filenames(j+12)
    load(filenames(j+12)) %sf 13 a 18
    SF_X_medio93(j,:) = 0.01*Phasor_mag_errmean';
    SF_X_std93(j,:) = 0.01*Phasor_mag_errstdev';
    SF_phi_medio93(j,:) = Phasor_ph_errmean';
    SF_phi_std93(j,:) = Phasor_ph_errstdev';
    SF_F_medio93(j,:) = 0.01*Freq_errmed';
    SF_F_std93(j,:) = 0.01*Freq_stddev';  
    
    filenames(j+18)
    load(filenames(j+18)) % sm 19 a 24
    SM_X_medio93(j,:) = 0.01*Phasor_mag_errmean';
    SM_X_std93(j,:) = 0.01*Phasor_mag_errstdev';
    SM_phi_medio93(j,:) = Phasor_ph_errmean';
    SM_phi_std93(j,:) = Phasor_ph_errstdev';
    SM_F_medio93(j,:) = 0.01*Freq_errmed';
    SM_F_std93(j,:) = 0.01*Freq_stddev';    
end

keepvars = {'SF_F_medio60', 'SF_F_medio93','SF_F_std60','SF_F_std93',...
    'SF_X_medio60', 'SF_X_medio93','SF_X_std60','SF_X_std93',...
'SF_phi_medio60', 'SF_phi_medio93','SF_phi_std60','SF_phi_std93',...
'SM_F_medio60', 'SM_F_medio93','SM_F_std60','SM_F_std93',...
    'SM_X_medio60', 'SM_X_medio93','SM_X_std60','SM_X_std93',...
'SM_phi_medio60', 'SM_phi_medio93','SM_phi_std60','SM_phi_std93',...
'filenames'};
clearvars('-except',keepvars{:})

ti = 1:9;
tau = ti/10;

fig_SF_F = figure(1); % fig salto fase - freq
subplot(2,1,1)
for j = 1:6
    eF = errorbar(tau,SF_F_medio60(j,:),SF_F_std60(j,:),SF_F_std60(j,:),'bx--');hold on;
end
grid on;
title ('Salto de Fase, SNR = 60 dB')
ylabel('Frequência estimada [Hz/Hz]')
subplot(212)
for j = 1:6
    eF = errorbar(tau,SF_F_medio93(j,:),SF_F_std93(j,:),SF_F_std93(j,:),'bx--');hold on;
end
grid on;
title ('Salto de Fase, SNR = 93 dB')
ylabel('Frequência estimada [Hz/Hz]')
xlabel('Localização relativa do salto')

fig_SF_X = figure(2)
subplot(2,1,1)
for j = 1:6
    eF = errorbar(tau,SF_X_medio60(j,:),SF_X_std60(j,:),SF_X_std60(j,:),'rx--');hold on;
end
title ('Salto de Fase, SNR = 60 dB')
ylabel('Magnitude estimada [V/V]')
grid on;
subplot(212)
for j = 1:6
    eF = errorbar(tau,SF_X_medio93(j,:),SF_X_std93(j,:),SF_X_std93(j,:),'rx--');hold on;
end
grid on;
title ('Salto de Fase, SNR = 93 dB')
ylabel('Magnitude estimada [V/V]')
xlabel('Localização relativa do salto')

fig_SF_phi = figure(3)
subplot(2,1,1)
for j = 1:6
    eF = errorbar(tau,SF_phi_medio60(j,:),SF_phi_std60(j,:),'mx--');hold on;
end
grid on;
title ('Salto de Fase, SNR = 60 dB')
ylabel('Fase estimada [graus]')
subplot(212)
for j = 1:6
    eF = errorbar(tau,SF_phi_medio93(j,:),SF_phi_std93(j,:),'mx--');hold on;
end
grid on;
title ('Salto de Fase, SNR = 93 dB')
ylabel('Fase estimada [graus]')
xlabel('Localização relativa do salto')

fig_SM_F = figure(4); % fig salto mg - freq
subplot(2,1,1)
for j = 1:6
    eF = errorbar(tau,SM_F_medio60(j,:),SM_F_std60(j,:),SM_F_std60(j,:),'bx--');hold on;
end
grid on;
title ('Salto de Magnitude, SNR = 60 dB')
ylabel('Frequência estimada [Hz/Hz]')
subplot(212)
for j = 1:6
    eF = errorbar(tau,SM_F_medio93(j,:),SM_F_std93(j,:),'bx--');hold on;
end
grid on;
title ('Salto de Magnitude, SNR = 93 dB')
ylabel('Frequência estimada [Hz/Hz]')
xlabel('Localização relativa do salto')

fig_SM_X = figure(5)
subplot(2,1,1)
for j = 1:6
    eF = errorbar(tau,SM_X_medio60(j,:),SM_X_std60(j,:),'rx--');hold on;
end
title ('Salto de Magnitude, SNR = 60 dB')
ylabel('Magnitude estimada [V/V]')
grid on;
subplot(212)
for j = 1:6
    eF = errorbar(tau,SM_X_medio93(j,:),SM_X_std93(j,:),'rx--');hold on;
end
grid on;
title ('Salto de Magnitude, SNR = 93 dB')
ylabel('Magnitude estimada [V/V]')
xlabel('Localização relativa do salto')


fig_SM_phi = figure(6)
subplot(2,1,1)
for j = 1:6
    eF = errorbar(tau,SM_phi_medio60(j,:),SM_phi_std60(j,:),'mx--');hold on;
end
grid on;
title ('Salto de Magnitude, SNR = 60 dB')
ylabel('Fase estimada [graus]')
subplot(212)
for j = 1:6
    eF = errorbar(tau,SM_phi_medio93(j,:),SM_phi_std93(j,:),'mx--');hold on;
end
grid on;
title ('Salto de Magnitude, SNR = 93 dB')
ylabel('Fase estimada [graus]')
xlabel('Localização relativa do salto')