%MC_HLM4
clear all; close all; clc;
MCiter = 1000;
SNR = 60; 
Fs = 5000;
runs = true;

%casos salto fase
hm = 0; ha = 10; phi0= 120; 
f=1;
% phase         X1  w    ph   x3 (KaS)
par_var_phase = [0  0 0 0]; % mantemos somente a incerteza de tau = 2dt
[filenames(f)] = MC_HLM4_par(MCiter,par_var_phase, SNR, Fs, hm, ha, phi0,runs)

%casos salto mag
hm = 0.1; ha = 0; phi0= 120;
% mag          x1  x2(KxS)  wf    ph  
par_var_mag = [0  0 0 0 ]; % parameter variation in percent related to nominal            
f = f+1;
[filenames(f)] = MC_HLM4_par(MCiter,par_var_mag, SNR, Fs, hm, ha, phi0,runs)

save('Cap4_LM\resultadosMC\filenames_par_tau.mat','filenames')
run('Cap4_LM\resultadosMC\analises_par_tau.m')