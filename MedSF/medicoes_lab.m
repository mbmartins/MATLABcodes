clear all; close all; clc;

% dados do ICHQP
load('ICHQP-DIP.mat')
%trecho com salto de fase/magnitude
Sig = Signal(1361:1361+499);

%trecho estavel com mag baixa
%Sig = Signal(1700:1700+499);

%trecho estavel com mag alta
%Sig = Signal(500:500+499);

plot(Sig);

%teriamos que detectar tau neste caso...
% primeira estimativa no olho
tau_n = 251;

[f1,f2,f_r,h_f] = MedSF(Sig, tau_n, Fs)

[f1p,f2p,f_rp,h_fp] = MedSF_PATV(Sig, tau_n, Fs)