%MC_HLM4
clear all; close all; clc;
MCiter = 1000;
SNR = 60; 
Fs = 5000;
runs = true;
f = 0;
F1 = 60;
Wr = 2*pi*F1;

%casos salto fase
hm = 0; ha = 10; phi0= 120; 
par_var_phase = [0 1 0 0]; % parameter variation in percent related to nominal
f = f+1;
[filenames(f), std_R2(f,:), std_w(f,:), std_Xmag(f,:), std_Xphi(f,:)] = MC_HLM4_par_w(MCiter,par_var_phase, SNR, Fs, hm, ha, phi0,runs)
ci_SF_w = max(std_w(f,:),[],2)/(2*Wr/sqrt(12))

%casos salto mag
hm = 0.1; ha = 0; phi0= 120;
% mag          x1  x2(KxS)  wf    ph  
par_var_phase = [0 0 1 0]; % parameter variation in percent related to nominal
f = f+1;
[filenames(f), std_R2(f,:), std_w(f,:), std_Xmag(f,:), std_Xphi(f,:)] = MC_HLM4_par_w(MCiter,par_var_phase, SNR, Fs, hm, ha, phi0,runs)         

ci_SM_w = max(std_w(f,:),[],2)

save('MC_ci_wconf')
% save('Cap4_LM\resultadosMC\filenames_par.mat','filenames')
% run('Cap4_LM\resultadosMC\analises_par.m')
