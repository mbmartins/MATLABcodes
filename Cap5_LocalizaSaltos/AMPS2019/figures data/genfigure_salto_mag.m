clear all; close all; clc;
%load('example_Signal.mat')

%gerar de novo
F0= 60;
F1 = 60;
SampleRate = 5000;
Ps = 0; %degrees
NCycles = 6;
tau1 = 0.7;
tau2 = 10;
SNR = 60;
ha = 0;
hx = 0.3;
KRS = 0;
Signal = SigGEN2(F0,F1,SampleRate,Ps,NCycles,tau1,tau2,SNR, ha, hx,KRS);
NSamples = length(Signal);

%ler já gerado
%load('sinal_exemplo_zero_cross.mat')

br = floor(0.05*NSamples);
z=hilbert(Signal');  % calculates the analytic signal associated with Signal
fi = unwrap(angle(z));
dt = 1/SampleRate;
T = NSamples*dt;
t = 0:dt:T-dt;


    df=gradient(fi);% Hilbert estimate of the instantaneous frequency of z
    df=abs(df-median(df(br:end-br))); %v3mod
    
    subplot(2,1,1)
    plot(df','.-'); ylabel("$|d[n]|$",'Interpreter','latex','FontSize',14)
    xlabel('$n$','Interpreter','latex','FontSize',14)
    grid on;
    subplot(2,1,2)
    plot(Signal','.-'); ylabel("$x[n]$",'Interpreter','latex','FontSize',14); 
        xlabel('$n$','Interpreter','latex','FontSize',14)
        grid on;