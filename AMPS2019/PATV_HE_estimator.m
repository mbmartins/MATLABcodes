function [tau_error,FE, dmax] = PATV_HE_estimator(SNR,KxS,KaS,Ps,tau1,SAG_cycles,lambda_a,lambda_theta)
% estimates frequency and tau_errors using NLHE estimator
%clear all; close all; clc;

% What differs from the HE_estimator???
% 1- Application of the PATV procedure to the signals:
% 1a) a_i[n], for the magnitude detector, zero order polinomial with steps
% The corresponding "clean" signal is 
% 1b) theta_i[n], for the phase detector, first order polinomial with steps

% 2 - Numeric differentiation of 1a) and 1b) provides detection signals
% thoroughly, either by the respective "u" signals or by the gradient of
% the respective "x" signals.


%nominal parameters
F0 = 60; F1 = 60; Fs = 4800;
%Ps = 90; 
NCycles = 6;
%SNR = 60;
%tau1 = 0.5;  % in [%] of the time window
%SAG_cycles = 1; %duration of SAG
%KaS = 0; %[degrees]
%KxS = -0.2; % [relative step]

th_a_i = 2.8e-3; th_fi = 4.8e-3; %thresholds for detection
%km = 2e7;
%kf = 100e9;

%OBS: the thresholds with PATV are not easy to set relative to the median,
%because the medians are usually too small.

br = 0.05; %percentage of samples to be ignored in detection 

NSamples = floor(NCycles*Fs/F0);
samples_cycle = Fs/F0;
br_mask = ones(1,NSamples); %mask to ignore samples
br_mask(1:floor(br*NSamples)) = 0; 
br_mask(floor((1-br)*NSamples):NSamples) = 0; 
tau2 = tau1+SAG_cycles*(Fs/F0)/NSamples; % time to end SAG in [%]

n = 1:NSamples;
th_a_i_v(1:NSamples) = th_a_i; %threshold for magnitude detection
th_fi_v(1:NSamples) = th_fi; %threshold for frequency detection

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

%loPATV params
% a = 1.0;
% % a = eps;    % L1 norm - reproduces above result
% phi_fi = @(x) lam_fi/a * log(1 + a*abs(x));
% phi_a_i = @(x) lam_a_i/a * log(1 + a*abs(x));
% %dphi = @(x) lam ./(1 + a*abs(x)) .* sign(x);
% wfun_fi = @(x) abs(x) .* (1 + a*abs(x)) / lam_fi;
% wfun_a_i = @(x) abs(x) .* (1 + a*abs(x)) / lam_a_i;

%PATV algorithm
%[x, p, cost, u, v] = patv_MM(y, d, lambda, Nit)
[x_a_i,p_a_i, cost_a_i, u_a_i, v_a_i] = patv_MM(a_i, d_a_i, lambda_a_i, Nit);
[x_theta_i, p_theta_i, cost_theta_i, u_theta_i, v_theta_i] = patv_MM(theta_i, d_theta_i, lambda_theta_i, Nit);

%loPATV algorithm
% [x_fi, p_theta_i, cost_theta_i, u, v] = patv_MM2(fi, d_fi,phi_fi,wfun_fi, Nit);
% [x_a_i,p_a_i, cost_a_i, u, v] = patv_MM2(a_i, d_a_i,phi_a_i,wfun_a_i, Nit);

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

% DETECTION SIGNALS
%grad_a_i = abs(gradient(x_a_i));
%grad_theta_i = abs(gradient(x_theta_i));
% OBS: in the PATV procedure, the output signal u is already calculated as
% a good approximation of the first order differences of x. 
grad_a_i = [abs(u_a_i - median(u_a_i)); 0];
grad_theta_i = [abs(u_theta_i - median(u_theta_i));0 ];

detector_a_i = br_mask'.*grad_a_i;
detector_theta_i = br_mask'.*grad_theta_i;

%limiar_mag=km*median(abs(detector_a_i));
%limiar_fase=kf*median(abs(detector_theta_i));
limiar_mag = th_a_i; 
limiar_fase = th_fi;


% DETECTION OF FIRST PEAK
[ga_i_max(1),imax_a_i(1)] = maxk(detector_a_i,1);
[gtheta_i_max(1),imax_theta_i(1)] = maxk(detector_theta_i,1); %detects first peak

%dmax = ga_i_max(1); % valor a ser retornado para tentar compensar o efeito sistemático
%dmax = median(abs(detector_a_i));
dmax = gtheta_i_max(1);

% subplot(2,1,2)
% plot(n,detector_a_i,'b');  ylabel('Detection signal')
%d2_a_i = abs(detector_a_i.*gradient(detector_a_i)); 
%hold on; plot(n,d2_a_i)

% DETECTION SIGNALS FOR THE SECOND PEAK
aind = floor(samples_cycle)/2; %half cycle number of samples 
detector_a_i2 = detector_a_i;
% guarantee indices within [1,NSamples]
n_ini = ((imax_a_i(1)-aind)>1)*(imax_a_i(1)-aind-1) + 1;
n_end = ((imax_a_i(1)+aind)<NSamples)*(imax_a_i(1)+aind-1) + 1;
%ignores half cycle centered in imax_a_i
detector_a_i2(n_ini:n_end) = 0;

detector_theta_i2 = detector_theta_i;
%ignores half cycle centered in imax_theta_i
n_ini = ((imax_theta_i(1)-aind)>1)*(imax_theta_i(1)-aind-1) + 1;
n_end = ((imax_theta_i(1)+aind)<NSamples)*(imax_theta_i(1)+aind-1) + 1;
detector_theta_i2(n_ini:n_end) = 0;

% DETECTION OF SECOND PEAK
[ga_i_max(2),imax_a_i(2)] = maxk(detector_a_i2,1);
[gtheta_i_max(2),imax_theta_i(2)] = maxk(detector_theta_i2,1);

% hold on; plot(n,detector_a_i,'red',n,th_a_i_v,'k--')
% legend('d1_{a_i}','d2_{a_i}', 'Location','southeast')

% figure(4)
% subplot(2,1,1)
% plot(n, fi,'.k', n, x_fi+p_theta_i, 'black') 
% title('Calculated TV component (PATV) - Frequency');
% legend('Data','Estimated signal', 'Location','southeast')
% subplot(2,1,2)
% plot(n,detector_fi,'b');  ylabel('Detection signal')

%ignore adjacent samples around the first peak
%detector_fi(imax_fi-br:imax_fi+br) = 0; 
% aind = imax_fi-floor((Fs/F0)/NSamples/2);
% while detector_theta_i(aind)>th_fi
%     detector_theta_i(aind:imax_fi) = 0;
%     aind = aind - 1;
% end
% aind = imax_fi+1;
% while detector_theta_i(aind)>th_fi
%     detector_theta_i(imax_fi:aind) = 0;
%     aind = aind + 1;
% end
% hold on; plot(n,detector_fi,'red',n,th_fi_v,'k--')
% legend('d1_{fi}','d2_{fi}', 'Location','southeast')


% error calculation
[tau_est_a_i, a_i_ind] = sort(imax_a_i); %error in samples
tau_est_a_i(ga_i_max(a_i_ind) < limiar_mag) = NaN; % check if estimates are valid
[tau_est_a_i,a_i_ind] = sort(tau_est_a_i(a_i_ind));
tau_error_a_i = (tau_est_a_i - [tau1 tau2]*NSamples); %error in [dt]

[tau_est_theta_i,theta_i_ind] = sort(imax_theta_i); %tau estimated in [samples]
tau_est_theta_i(gtheta_i_max(theta_i_ind) < limiar_fase) = NaN; % check if estimates are valid
[tau_est_theta_i,theta_i_ind] = sort(tau_est_theta_i(theta_i_ind));
tau_error_theta_i = (tau_est_theta_i - [tau1 tau2]*NSamples); %error in [dt]

%criterio de escolha mag x freq com denoising: usar a maior diferença para 
% o limiar ??
%crit = ga_i_max(a_i_ind) - th_a_i > gtheta_i_max(theta_i_ind) - th_fi;
%% OU
%% usar a razão, para manter o paralelismo

crit = ga_i_max(a_i_ind)/limiar_mag > gtheta_i_max(theta_i_ind)/limiar_fase;
tau_est(crit) = tau_est_a_i(crit);
tau_est(~crit) = tau_est_theta_i(~crit);
tau_error = tau_est -2 - [tau1 tau2]*NSamples;


if abs(tau_error(1))>2
    %false positive
    ratio_a = ga_i_max(a_i_ind)/limiar_mag;
    ratio_f = gtheta_i_max(theta_i_ind)/limiar_fase;
    if ratio_a > ratio_f
        tau_est;
    end
end
if isnan(tau_error(1))
    tau_est;
end



%calcular frequencia pela inclinação da fase
%P = polyfit(n,p_theta_i',1);
%f_est = P(1)*Fs/(2*pi);

%calcular a freq pela mediana da freq istantanea
f_est = median(diff(p_theta_i))*Fs/(2*pi);

FE = (f_est - F1); %Hz


tau_est = (tau_est)/Fs;


%problemas: 
% 1 - indices excedendo as bordas (na hora de excluir do sinal de
% detecção) -> implementar um teste para isso não ocorrer
% por exemplo, baseado no th_a_i - ok, mas usando meio ciclo em torno do
% pico funcionou melhor

% 2 - picos ocorrendo perto das bordas (falsos ou verdadeiros) quando f1
% está fora da nominal
% -> implementar exclusão próximo das bordas (2 a 5%) ??
% -> picos verdadeiros terão que ser excluídos -> limitação do método
% talvez identificar de alguma forma o número de degraus pelo PATV
% o PATV "filtra" as não idealidades de Hilbert próximo a borda

% 3 - redução do desempenho da detecção quando os degraus estão próximos da
% borda -> limitação do método -> determinar os limites

% 4 - quando só existe tau2, fazendo tau1<=0 -> ver como o sort lida com
% NaN -> na verdade há um problema no cálculo do erro...

% 5 - possibilidade de usar o LoPATV para melhoria do desempenho, se
% necessário - não funcionou bem...
