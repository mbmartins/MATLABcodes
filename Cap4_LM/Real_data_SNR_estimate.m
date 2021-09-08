% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

D = open('Resultados_30_05.mat')
D2 = open('Steps_complete.mat')

dt = 2.0e-4;
SampleRate = 1/dt; %Hz

Signal = [
 D.DadosS1.Cm60_5000';
 D2.Dadost5.MS_p_0';
 ];

x = Signal(1,1:2000);
snr(x,SampleRate)

