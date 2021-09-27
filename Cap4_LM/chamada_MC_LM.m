%MC_HLM4
clear all; close all; clc;
MCiter = 1000;
SNR = [50]; 
Fs = 5000;
runs = true;
f = 0;
for s = 1
%casos salto fase
hm = 0; ha = 10; phi0= 0; f = f+1;
[filenames(f)] = MC_HLM4(MCiter, SNR(s), Fs, hm, ha, phi0,runs)
hm = 0; ha = 10; phi0= 120; f = f+1;
[filenames(f)] = MC_HLM4(MCiter, SNR(s), Fs, hm, ha, phi0,runs)
hm = 0; ha = 10; phi0= -120; f = f+1;
[filenames(f)] = MC_HLM4(MCiter, SNR(s), Fs, hm, ha, phi0,runs)
hm = 0; ha = -10; phi0= 0; f = f+1;
[filenames(f)] = MC_HLM4(MCiter, SNR(s), Fs, hm, ha, phi0,runs)
hm = 0; ha = -10; phi0= 120; f = f+1;
[filenames(f)] = MC_HLM4(MCiter, SNR(s), Fs, hm, ha, phi0,runs)
hm = 0; ha = -10; phi0= -120; f = f+1;
[filenames(f)] = MC_HLM4(MCiter, SNR(s), Fs, hm, ha, phi0,runs)
% num total de 6, são 6x9000 = 54000 rodadas para salto fase

ha = 0;
%casos salto magnitude
hm = 0.1;  phi0= 0;f = f+1;
[filenames(f)] = MC_HLM4(MCiter, SNR(s), Fs, hm, ha, phi0,runs)
hm = 0.1; ha = 0; phi0= 120; f = f+1;
[filenames(f)] = MC_HLM4(MCiter, SNR(s), Fs, hm, ha, phi0,runs)
hm = 0.1; phi0= -120; f = f+1;
[filenames(f)] = MC_HLM4(MCiter, SNR(s), Fs, hm, ha, phi0,runs)
hm = -0.1; phi0= 0;f = f+1;
[filenames(f)] = MC_HLM4(MCiter, SNR(s), Fs, hm, ha, phi0,runs)
hm = -0.1; phi0= 120;f = f+1;
[filenames(f)] = MC_HLM4(MCiter, SNR(s), Fs, hm, ha, phi0,runs)
hm = -0.1; phi0= -120; f = f+1;
[filenames(f)] = MC_HLM4(MCiter, SNR(s), Fs, hm, ha, phi0,runs)
end

save('MATLABcodes\Cap4_LM\resultadosMC\filenames.mat','filenames')
% run('Cap4_LM\resultadosMC\analises.m')
run('MATLABcodes\Cap4_LM\resultadosMC\analises50.m')