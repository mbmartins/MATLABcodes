clear all; close all; clc;

F0 = 60;
Fs = 4800;
Ps = 90;
SNR = 60;

Signal = SigGen(F0,Fs,Ps,NCycles,tau1,tau2,SNR);