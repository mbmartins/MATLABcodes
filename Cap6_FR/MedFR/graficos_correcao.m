clear all; close all; clc;

load('analise_correcao')
subplot(211)
p1 = plot(phi_0_mean(1:end-2),Cf(1:end-2),'b.-')%,phi_0_mean(1:end-2),Cf(1:end-2))
xlim([-10 180])
xlb = xlabel({'$\hat{\phi}_0$','a)'}); 
xlb.Interpreter = 'latex'
ylb = ylabel('C_f [Hz/Hz]'); grid on
subplot(212)
p2 = plot(phi_0_mean(1:end-2)+5,dCf(1:end-2),'b.-'); grid on
%p2.axis = [-20 180 -1.5e-5 1.5e-5];
xlim([-10 180]); ylim([-1.5e-5 1.5e-5])
xlb = xlabel({'$\hat{\phi}_0$','b)'}); 
xlb.Interpreter = 'latex'
ylabel('\Delta C_f / \Delta \phi_0')


% figure(2)
% plot(Ps(1:end-2),Ephi(1:end-2));

