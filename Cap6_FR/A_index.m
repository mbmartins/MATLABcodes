% Estimadores de frequencia

run('teste')
save('teste')

%salto fase
run('estimadoresFR_salto_fase_phi0'); save('salto fase\phi0')
run('estimadoresFR_salto_fase_tau_n'); save('salto fase\tau_n')
run('estimadoresFR_salto_fase_T'); save('salto fase\T')
run('estimadoresFR_salto_fase_Fs'); save('salto fase\Fs')
run('estimadoresFR_salto_fase_F1'); save('salto fase\F1')
run('estimadoresFR_salto_fase_kf'); save('salto fase\ka')

%salto mag
run('estimadoresFR_salto_mag_phi0'); save('salto mag\phi0')
run('estimadoresFR_salto_mag_tau_n'); save('salto mag\tau_n')
run('estimadoresFR_salto_mag_T'); save('salto mag\T')
run('estimadoresFR_salto_mag_Fs'); save('salto mag\Fs')
run('estimadoresFR_salto_mag_F1'); save('salto mag\F1')
run('estimadoresFR_salto_mag_kf'); save('salto mag\kx')
