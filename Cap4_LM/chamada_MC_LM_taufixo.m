%MC_HLM4_taufixo
clear all; close all; clc;
MCiter = 1000;
SNR = [93.5]; 
Fs = 5000;
runs = true;
f = 0;
ti = 5;
s = 1;
%casos salto fase
hm = 0; ha = 10; phi0= 360; f = f+1;
[filenames(f)] = MC_HLM4_taufixo(MCiter, SNR(s), Fs, hm, ha, phi0,runs,ti)
hm = 0; ha = 10; phi0= 120; f = f+1;
[filenames(f)] = MC_HLM4_taufixo(MCiter, SNR(s), Fs, hm, ha, phi0,runs,ti)
hm = 0; ha = 10; phi0= -120; f = f+1;
[filenames(f)] = MC_HLM4_taufixo(MCiter, SNR(s), Fs, hm, ha, phi0,runs,ti)
hm = 0; ha = -10; phi0= 360; f = f+1;
[filenames(f)] = MC_HLM4_taufixo(MCiter, SNR(s), Fs, hm, ha, phi0,runs,ti)
hm = 0; ha = -10; phi0= 120; f = f+1;
[filenames(f)] = MC_HLM4_taufixo(MCiter, SNR(s), Fs, hm, ha, phi0,runs,ti)
hm = 0; ha = -10; phi0= -120; f = f+1;
[filenames(f)] = MC_HLM4_taufixo(MCiter, SNR(s), Fs, hm, ha, phi0,runs,ti)
% num total de 6, são 6x9000 = 54000 rodadas para salto fase

ha = 0;
%casos salto magnitude
hm = 0.1;  phi0= 360;f = f+1;
[filenames(f)] = MC_HLM4_taufixo(MCiter, SNR(s), Fs, hm, ha, phi0,runs,ti)
hm = 0.1; ha = 0; phi0= 120; f = f+1;
[filenames(f)] = MC_HLM4_taufixo(MCiter, SNR(s), Fs, hm, ha, phi0,runs,ti)
hm = 0.1; phi0= -120; f = f+1;
[filenames(f)] = MC_HLM4_taufixo(MCiter, SNR(s), Fs, hm, ha, phi0,runs,ti)
hm = -0.1; phi0= 360;f = f+1;
[filenames(f)] = MC_HLM4_taufixo(MCiter, SNR(s), Fs, hm, ha, phi0,runs,ti)
hm = -0.1; phi0= 120;f = f+1;
[filenames(f)] = MC_HLM4_taufixo(MCiter, SNR(s), Fs, hm, ha, phi0,runs,ti)
hm = -0.1; phi0= -120; f = f+1;
[filenames(f)] = MC_HLM4_taufixo(MCiter, SNR(s), Fs, hm, ha, phi0,runs,ti)


save('Cap4_LM\resultadosMC\filenames.mat','filenames')
run('Cap4_LM\resultadosMC\analises_taufixo.m')