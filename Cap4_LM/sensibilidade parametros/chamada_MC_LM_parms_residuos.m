%MC_HLM4
clear all; close all; clc;
MCiter = 1000;
SNR = 60; 
Fs = 5000;
runs = true;
f = 0;

%casos salto fase
hm = 0; ha = 10; phi0= 120; f = f+1;
% phase         X1  w    ph   x3 (KaS)
par_var_phase = [1  0 0 0]; % parameter variation in percent related to nominal
[filenames(f), std_R2(f,:), std_w(f,:), std_Xmag(f,:), std_Xphi(f,:)] = MC_HLM4_par_Xm(MCiter,par_var_phase, SNR, Fs, hm, ha, phi0,runs)
par_var_phase = [0 1 0 0]; % parameter variation in percent related to nominal
f = f+1;
[filenames(f), std_R2(f,:), std_w(f,:), std_Xmag(f,:), std_Xphi(f,:)] = MC_HLM4_par_w(MCiter,par_var_phase, SNR, Fs, hm, ha, phi0,runs)
par_var_phase = [0 0 1 0]; % parameter variation in percent related to nominal
f = f+1;
[filenames(f), std_R2(f,:), std_w(f,:), std_Xmag(f,:), std_Xphi(f,:)] = MC_HLM4_par_phi0(MCiter,par_var_phase, SNR, Fs, hm, ha, phi0,runs)
par_var_phase = [0 0 0 1]; % parameter variation in percent related to nominal
f = f+1;
[filenames(f), std_R2(f,:), std_w(f,:), std_Xmag(f,:), std_Xphi(f,:)] = MC_HLM4_par_hx(MCiter,par_var_phase, SNR, Fs, hm, ha, phi0,runs)


%casos salto mag
hm = 0.1; ha = 0; phi0= 120;
% mag          x1  x2(KxS)  wf    ph  
par_var_mag = [1  0 0 0 ]; % parameter variation in percent related to nominal            
f = f+1;
[filenames(f), std_R2(f,:), std_w(f,:), std_Xmag(f,:), std_Xphi(f,:)] = MC_HLM4_par_Xm(MCiter,par_var_phase, SNR, Fs, hm, ha, phi0,runs)
par_var_mag = [0 1 0 0 ]; % parameter variation in percent related to nominal            
f = f+1;
[filenames(f), std_R2(f,:), std_w(f,:), std_Xmag(f,:), std_Xphi(f,:)] = MC_HLM4_par_hx(MCiter,par_var_phase, SNR, Fs, hm, ha, phi0,runs)
par_var_mag = [0 0 1 0 ]; % parameter variation in percent related to nominal            
f = f+1;
[filenames(f), std_R2(f,:), std_w(f,:), std_Xmag(f,:), std_Xphi(f,:)] = MC_HLM4_par_w(MCiter,par_var_phase, SNR, Fs, hm, ha, phi0,runs)
par_var_mag = [0  0 0 1]; % parameter variation in percent related to nominal            
f = f+1;
[filenames(f), std_R2(f,:), std_w(f,:), std_Xmag(f,:), std_Xphi(f,:)] = MC_HLM4_par_phi0(MCiter,par_var_phase, SNR, Fs, hm, ha, phi0,runs)

Wr = 2*pi*60; Xr = 1; Xphir = 0.01*120;
% corrigir para valores relativos 
ci_SF_R2 = max(std_R2(1:4,:),[],2)
ci_SF_w = max(std_w(1:4,:),[],2)/(2*Wr/sqrt(12))
ci_SF_Xmag = max(std_Xmag(1:4,:),[],2)/(2*Xr/sqrt(12))
ci_SF_Xphi = max(std_Xphi(1:4,:),[],2)/(2*Xphir/sqrt(12))

ci_SM_R2 = max(std_R2(5:8,:),[],2)
ci_SM_w = max(std_w(5:8,:),[],2)/(2*Wr/sqrt(12))
ci_SM_Xmag = max(std_Xmag(5:8,:),[],2)/(2*Xr/sqrt(12))
ci_SM_Xphi = max(std_Xphi(5:8,:),[],2)/(2*Xphir/sqrt(12))

save('Cap4_LM\sensibilidade parametros\MC_ci100')
% save('Cap4_LM\resultadosMC\filenames_par.mat','filenames')
% run('Cap4_LM\resultadosMC\analises_par.m')
