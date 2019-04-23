clear all; close all; clc;

%nominal parameters
F0 = 60; F1 = 60;
Fs = 4800;
Ps = 90;
NCycles = 6;
NSamples = floor(NCycles*Fs/F0);
n = 1:NSamples;

%parameters to set
SNR = 90;
tau1 = 0.5;  % in [%] of the time window
tau2 = tau1+(Fs/F0)/NSamples;  %+1 cycle
KaS = 10; %[degrees]
KxS = -0.10; % [relative step]
th_gmi = 1e-3; th_fi = 1e-3;
th_gmi_v(1:NSamples) = th_gmi; %threshold for magnitude detection
th_fi_v(1:NSamples) = th_fi; %threshold for frequency detection

Signal = SigGEN(F0,F1,Fs,Ps,NCycles,tau1,tau2,SNR,KaS, KxS);
%SigGEN(F0,F1,SampleRate,Ps,NCycles,tau1,tau2,SNR, KaS, KxS)

%calculate instantaneous freq and magnitude by Hilbert Transform
z=hilbert(Signal');  % calculates the analytic signal associated with Signal
fi = unwrap(angle(z)); gmi = abs(z);

figure(1)

subplot(3,1,1); plot(Signal); title('Signal')
subplot(3,1,2); plot(fi); title('Instantaneous frequency')
subplot(3,1,3); plot(gmi); title('Instantaneous magnitude')


%fazer o denoising por PATV
%% Perform PATV filtering
% PATV: Least square polynomial approximation + total variation denoising

% parameters
d_fi = 1;                          % d : degree of approximation polynomial
d_gmi = 0;
lam_fi = 1;
lam_gmi = 0.5;   % lambda : regularization parameter
Nit = 10;                      % Nit : number of iterations

[x_fi, p_fi, cost_fi, u, v] = patv_MM(fi, d_fi, lam_fi, Nit);
[x_gmi,p_gmi, cost_gmi, u, v] = patv_MM(gmi, d_gmi, lam_gmi, Nit);

% display cost function history
figure(2)
subplot(2,1,1)
semilogy(1:Nit, cost_fi)
title('PATV algorithm - Cost function history');
xlabel('Iteration')
ylabel('Cost function value')
subplot(2,1,2)
semilogy(1:Nit, cost_gmi)
title('PATV algorithm - Cost function history');
xlabel('Iteration')
ylabel('Cost function value')

%display TV component
figure(3)
clf
% subplot(3,1,1)
% plot(n, fi,'.k', n, x_fi+p_fi, 'black') 
% title('Calculated TV component (PATV) - Frequency');
% legend('Data','Estimated signal', 'Location','southeast')
subplot(2,1,1)
plot(n, gmi,'.k', n, x_gmi+p_gmi, 'black') 
title('Calculated TV component (PATV) - Magnitude');
legend('Data','Estimated signal', 'Location','southeast')
detector_gmi = abs(gradient(x_gmi));
subplot(2,1,2)
plot(n,detector_gmi,'b');  ylabel('Detection signal')
%d2_gmi = abs(detector_gmi.*gradient(detector_gmi)); hold on; plot(n,d2_gmi)

%%%%%%%%%%%  estimate taus

%detector_gmi
[gmi_max,imax_gmi(1)] = maxk(detector_gmi,1);
%[gmi_min,imin] = min(detector_gmi);
% ignores adjacent samples
br = floor(0.02*NSamples);
detector_gmi(imax_gmi-br:imax_gmi+br) = 0;
hold on; plot(n,detector_gmi,'red',n,th_gmi_v,'k--')
legend('d1_{gmi}','d2_{gmi}', 'Location','southeast')
[gmi_max2,imax_gmi(2)] = maxk(detector_gmi,1);
tau_est_gmi = sort(imax_gmi); %error in samples
tau1_error_gmi = (tau_est_gmi(1) - tau1*NSamples) %error in [dt]
tau2_error_gmi = (tau_est_gmi(2) - tau2*NSamples) %error in [dt]

%detector_fi
detector_fi = abs(gradient(x_fi));
figure(4)
subplot(2,1,1)
plot(n, fi,'.k', n, x_fi+p_fi, 'black') 
title('Calculated TV component (PATV) - Frequency');
legend('Data','Estimated signal', 'Location','southeast')
subplot(2,1,2)
plot(n,detector_fi,'b');  ylabel('Detection signal')

[fi_max,imax_fi(1)] = maxk(detector_fi,1); %detects first tau
detector_fi(imax_fi-br:imax_fi+br) = 0; %ignores adjacent samples
hold on; plot(n,detector_fi,'red',n,th_fi_v,'k--')
legend('d1_{fi}','d2_{fi}', 'Location','southeast')
[fi_max2,imax_fi(2)] = maxk(detector_fi,1);
tau_est_fi = sort(imax_fi); %error in samples
tau1_error_fi = (tau_est_fi(1) - tau1*NSamples) %error in [dt]
tau2_error_fi = (tau_est_fi(2) - tau2*NSamples) %error in [dt]

% 1 - check if estimates are valid
if (tau_est_fi(1) <= th_fi)
    tau_est_fi = NaN;

% individualmente ou em grupo??
crit_fi = fi_max/th_fi;
crit_gmi = gmi_max/th_gmi;

%calcular frequencia pela inclinação da fase

