function [FE1, FE2, FE3] = PATV_FE_estimator(SNR,KxS,KaS,Ps,tau1,SAG_cycles,lambda_a,lambda_theta)
% estimates frequency NLHE estimator
%clear all; close all; clc;

%nominal parameters
F0 = 60; F1 = 60; Fs = 4800;
%Ps = 90; 
NCycles = 6;
%SNR = 60;
%tau1 = 0.5;  % in [%] of the time window
%SAG_cycles = 1; %duration of SAG
%KaS = 0; %[degrees]
%KxS = -0.2; % [relative step]

br = 0.05; %percentage of samples to be ignored in detection 

NSamples = floor(NCycles*Fs/F0);
samples_cycle = Fs/F0;
br_mask = ones(1,NSamples); %mask to ignore samples
br_mask(1:floor(br*NSamples)) = 0; 
br_mask(floor((1-br)*NSamples):NSamples) = 0; 
tau2 = tau1+SAG_cycles*(Fs/F0)/NSamples; % time to end SAG in [%]

n = 1:NSamples;

Signal = SigGEN(F0,F1,Fs,Ps,NCycles,tau1,tau2,SNR,KaS, KxS);
%SigGEN(F0,F1,SampleRate,Ps,NCycles,tau1,tau2,SNR, KaS, KxS)

%calculate instantaneous freq and magnitude by Hilbert Transform
z=hilbert(Signal');  % calculates the analytic signal associated with Signal
theta_i = unwrap(angle(z)); a_i = abs(z);

% figure(1)
% plot(Signal,'.k'); ylabel('Sampled signal x[n] [V]'); xlabel('Samples')
% axis([1 480 -1 1]); grid on;

% subplot(3,1,1); plot(Signal); title('Signal'); grid on
% subplot(3,1,2); plot(fi); title('Instantaneous phase'); xlabel('Samples');ylabel('\theta_i[n]')
% subplot(3,1,3); plot(a_i); title('Instantaneous magnitude'); xlabel('Samples');ylabel('\a_i[n]')

%fazer o denoising por PATV
%% Perform PATV filtering
% PATV: Least square polynomial approximation + total variation denoising

% parameters
d_theta_i = 1;                          % d : degree of approximation polynomial
d_a_i = 0;
lambda_a_i = lambda_a;   % lambda : regularization parameter
lambda_theta_i = lambda_theta;
Nit = 10;                      % Nit : number of iterations

%PATV algorithm
%[x, p, cost, u, v] = patv_MM(y, d, lambda, Nit)

% FE1 - patv em theta_i
[x_theta_i, p_theta_i, cost_theta_i, u_theta_i, v_theta_i] = patv_MM(theta_i, d_theta_i, lambda_theta_i, Nit);
% FE2 - patv em f_i
f_i = gradient(theta_i)*Fs/(2*pi); d=0;
[x_theta_i2, p_f_i, cost_theta_i2, u_theta_i2, v_theta_i2] = patv_MM(f_i, d, lambda_theta_i, Nit);

% display cost function history
% figure(2)
% subplot(2,1,1)
% semilogy(1:Nit, cost_theta_i)
% title('PATV algorithm - Cost function history');
% xlabel('Iteration')
% ylabel('Cost function value')
% subplot(2,1,2)
% semilogy(1:Nit, cost_a_i)
% title('PATV algorithm - Cost function history');
% xlabel('Iteration')
% ylabel('Cost function value')

%display TV component
% figure(3)
% clf
% % subplot(3,1,1)
% % plot(n, fi,'.k', n, x_fi+p_theta_i, 'black') 
% % title('Calculated TV component (PATV) - Frequency');
% % legend('Data','Estimated signal', 'Location','southeast')
% subplot(2,1,1)
% plot(n, a_i,'.k', n, x_a_i+p_a_i, 'black') 
% title('Calculated TV component (PATV) - Magnitude');
% legend('Data','Estimated signal', 'Location','southeast')

%calcular a freq pela mediana da freq instantanea
f_i1 = gradient(p_theta_i)*Fs/(2*pi);
f_est1 = median(f_i1);

% mediana de f já filtrado pelo PATV
f_est2 = median(p_f_i);

%other PATV for frequency estimation
lambda3 = 1.;
[x_theta_i, p3_theta_i, cost_theta_i, u_theta_i, v_theta_i] = patv_MM(theta_i, d_theta_i, lambda3, Nit);
f_i3 = gradient(p3_theta_i)*Fs/(2*pi);
f_est3 = median(f_i3);




phi_0_est = p_theta_i(1)*180/pi

%fazer função para correção sistemática
%freq_corr = 0.00298276140145513; %phi0 = 0
% freq_corr = 0.00296409387030716; %phi0 = 15
% freq_corr = 0.00296598890104022; %phi0 = 30
% freq_corr = 0.00306010061155422; %phi0 = 45
% freq_corr = 0.00323290947446844; %phi0 = 60
% freq_corr = 0.00338467357633335; %phi0 = 75
% freq_corr = 0.00342652604283979; %phi0 = 90
freq_corr = 0;

FE1 = (f_est1 - F1 - freq_corr*F1); %Hz

FE2 = (f_est2 - F1 - freq_corr*F1); %Hz

FE3 = (f_est3 - F1 - freq_corr*F1); %Hz

