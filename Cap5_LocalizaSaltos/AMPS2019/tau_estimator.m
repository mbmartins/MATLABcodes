function [tau1,tau2,f_est] = tau_estimator(Signal,F0,Fs)
samples_cycle = Fs/F0;
NSamples = length(Signal);
n = 1:NSamples;
%possible parameters for function, now hardcoded
th_gmi = 1e-3; th_fi = 1e-4; %thresholds for detection
br = 0.02; %percentage of samples to be ignored in detection

%calculate instantaneous freq and magnitude by Hilbert Transform
z=hilbert(Signal');  % calculates the analytic signal associated with Signal
fi = unwrap(angle(z)); gmi = abs(z);

figure(1)
subplot(3,1,1); plot(Signal); title('Signal'); grid on
subplot(3,1,2); plot(fi); title('Instantaneous frequency')
subplot(3,1,3); plot(gmi); title('Instantaneous magnitude')

% PATV: Least square polynomial approximation + total variation denoising
% parameters
d_fi = 1;                          % d : degree of approximation polynomial
d_gmi = 0;
lam_fi = 1;
lam_gmi = 0.5;   % lambda : regularization parameter
Nit = 10;                      % Nit : number of iterations

%PATV algorithm
[x_fi, p_fi, cost_fi, u, v] = patv_MM(fi, d_fi, lam_fi, Nit);
[x_gmi,p_gmi, cost_gmi, u, v] = patv_MM(gmi, d_gmi, lam_gmi, Nit);

%%% Hybrid detectors
br_mask = ones(1,NSamples); %mask to ignore samples
br_mask(1:floor(br*NSamples)) = 0; 
br_mask(floor((1-br)*NSamples):NSamples) = 0; 
%display TV component
figure(3)
th_gmi_v(1:NSamples) = th_gmi; %threshold for magnitude detection
th_fi_v(1:NSamples) = th_fi; %threshold for frequency detection
subplot(3,1,1)
plot(n,Signal,'.k'); grid on; xlabel('Samples'); ylabel('x[n] [V]');
title('Sampled signal');

subplot(3,1,2)
plot(n, gmi,'.k', n, x_gmi+p_gmi, 'black') 
title('Instantaneous Magnitude');
legend('a_i[n]','w[n]', 'Location','southeast')
detector_gmi = br_mask'.*abs(gradient(x_gmi));
xlabel('Samples'); ylabel('w[n]')

subplot(3,1,3)
plot(n,detector_gmi); xlabel('Samples'); 

%%%%%%%%%%%  estimate taus

%detector_gmi
[gmi_max(1),imax_gmi(1)] = maxk(detector_gmi,1);
%[gmi_min,imin] = min(detector_gmi);
% ignores adjacent samples
%br = floor(0.02*NSamples);
aind = floor(samples_cycle/3); %ignores half cycle centered in imax_gmi
detector_gmi(imax_gmi(1)-aind:imax_gmi(1)+aind) = 0;
%aind = imax_gmi-1;
% while detector_gmi(aind)>th_gmi
%     detector_gmi(aind:imax_gmi) = 0;
%     aind = aind - 1;
% end
% aind = imax_gmi+1;
% while detector_gmi(aind)>th_gmi
%     detector_gmi(imax_gmi:aind) = 0;
%     aind = aind + 1;
% end
hold on; plot(n,detector_gmi,n,th_gmi_v,'k--')
title('Detection signals')
legend('d1_{ai}','d2_{ai}', 'Threshold', 'Location','southeast')
xlabel('Samples'); ylabel('d_{mNL}[n]')
[gmi_max(2),imax_gmi(2)] = maxk(detector_gmi,1);

%detector_fi
detector_fi = br_mask'.*abs(gradient(x_fi));
figure(4)
subplot(2,1,1)
plot(n, fi,'.k', n, x_fi+p_fi, 'black') 
title('Calculated TV component (PATV) - Frequency');
legend('Data','Estimated signal', 'Location','southeast')
subplot(2,1,2)
plot(n,detector_fi,'b');  ylabel('Detection signal')

[fi_max(1),imax_fi(1)] = maxk(detector_fi,1); %detects first peak

%ignore adjacent samples
%detector_fi(imax_fi-br:imax_fi+br) = 0; 
aind = imax_fi-0.5*(floor((Fs/F0)/NSamples));
while detector_fi(aind)>th_fi
    detector_fi(aind:imax_fi) = 0;
    aind = aind - 1;
end
aind = imax_fi+1;
while detector_fi(aind)>th_fi
    detector_fi(imax_fi:aind) = 0;
    aind = aind + 1;
end
hold on; plot(n,detector_fi,'red',n,th_fi_v,'k--')
legend('d1_{fi}','d2_{fi}', 'Location','southeast')
[fi_max(2),imax_fi(2)] = maxk(detector_fi,1); %detects second peak

% tau calculation
[tau_est_gmi, gmi_ind] = sort(imax_gmi); %error in samples
tau_est_gmi(gmi_max < th_gmi) = NaN; % check if estimates are valid
[tau_est_fi,fi_ind] = sort(imax_fi); %tau estimated in [samples]
tau_est_fi(fi_max < th_fi) = NaN; % check if estimates are valid

%criterio de escolha mag x freq com denoising: usar a maior diferença para 
% o limiar
crit = gmi_max(gmi_ind) - th_gmi > fi_max(fi_ind) - th_fi;
tau_est(crit) = tau_est_gmi(crit)
tau_est(~crit) = tau_est_fi(~crit)
%tau_error = tau_est - [tau1 tau2]*NSamples

%calcular frequencia pela inclinação da fase
P = polyfit(n,p_fi',1);
f_est = P(1)*Fs/(2*pi);
%FE = (f_est - F1)*100/F1

tau1 = tau_est(1);
tau2 = tau_est(2);