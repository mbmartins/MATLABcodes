function [f_r, phi_0] = MedFR_PATV(Signal,lambda,Fs)
% estimates frequency 
 
%nominal parameters
%F0 = 60; F1 = 60; Fs = 4800;
%Ps = 90; 
%NCycles = 6;

%br = 0.05; %percentage of samples to be ignored in detection 
%NSamples = floor(NCycles*Fs/F0);
%samples_cycle = Fs/F0;
%br_mask(1:floor(br*NSamples)) = 0; 
%br_mask(floor((1-br)*NSamples):NSamples) = 0; 
%tau2 = tau1+SAG_cycles*(Fs/F0)/NSamples; % time to end SAG in [%]

%n = 1:NSamples;

%calculate instantaneous freq and magnitude by Hilbert Transform
z=hilbert(Signal');  % calculates the analytic signal associated with Signal
theta_i = unwrap(angle(z)); %a_i = abs(z);

% parameters
d_theta_i = 1;                          % d : degree of approximation polynomial
Nit = 20;                      % Nit : number of iterations

%PATV algorithm
[~, p_theta_i, ~, ~, ~] = patv_MM(theta_i, d_theta_i, lambda, Nit);
phi_0 = p_theta_i(1)*180/pi; % in degrees
%grad_theta_i = [abs(u_theta_i);0 ];

%calcular a freq pela mediana da freq instantanea
f_i = gradient(p_theta_i)*Fs/(2*pi);

f_r = median(f_i);

% %fazer função para correção sistemática
% %freq_corr = 0.00298276140145513; %phi0 = 0
% % freq_corr = 0.00296409387030716; %phi0 = 15
% % freq_corr = 0.00296598890104022; %phi0 = 30
% % freq_corr = 0.00306010061155422; %phi0 = 45
% % freq_corr = 0.00323290947446844; %phi0 = 60
% % freq_corr = 0.00338467357633335; %phi0 = 75
% % freq_corr = 0.00342652604283979; %phi0 = 90
% freq_corr = 0;

