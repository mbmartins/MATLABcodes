% Estimation of Frequency of lab signals
clear all; close all; clc;

%signal from memory
load('60Hz_DMM.mat');
NSamples = size(Signal,1);
F1 = 60; %nominal frequency [Hz]

%Hilbert filter
br = floor(0.05*NSamples); % 5% of NSamples are taken off at the beggining and end
z=hilbert(Signal);  % calculates the analytic signal associated with Signal
Psi = unwrap(angle(z));

f_hi=unwrap(angle(z(2:end,:).*conj(z(1:end-1,:))));  % Hilbert estimate of the instantaneous frequency of Signal
f=f_hi-median(f_hi(br:end-br));

plot(t,Signal)

%estimation of angular frequency (w = 2*pi*f) using hilbert
% w is the slope of the linear curve

%1st degree model matrix
%calculating 2freqs
X1 = [t ones(size(t))]; y = Psi;
beta = (X1'*X1)\X1'*y;
fest_1 = beta(1)/(2*pi)
ferr_1 = 100*(fest_1 - F1)/F1
res_1 = (X1*beta - y);

figure
subplot(3,1,1); plot(t,Psi); title('\Psi')
subplot(3,1,2); plot(t,res_1,'.'); title('residuals = (X*\beta - y)')
subplot(3,1,3); plot(t(2:end),abs(f_hi)); title('Instantaneous frequency')


