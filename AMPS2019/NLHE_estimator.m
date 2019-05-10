function [FE,tau_error] = HE(SNR,KxS,KaS,Ps,tau1,SAG_cycles)
% estimates frequency and number of tau_error > 2dt using HE estimator
%clear all; close all; clc;

%nominal parameters
F0 = 60; F1 = 60; Fs = 4800;
Ps = 90; 
NCycles = 6;
%SNR = 60;
%tau1 = 0.5;  % in [%] of the time window
%SAG_cycles = 1; %duration of SAG
%KaS = 0; %[degrees]
%KxS = -0.2; % [relative step]

th_gmi = 1e-3; th_fi = 1e-4; %thresholds for detection
br = 0.05; %percentage of samples to be ignored in detection 

NSamples = floor(NCycles*Fs/F0);
samples_cycle = Fs/F0;
br_mask = ones(1,NSamples); %mask to ignore samples
br_mask(1:floor(br*NSamples)) = 0; 
br_mask(floor((1-br)*NSamples):NSamples) = 0; 
tau2 = tau1+SAG_cycles*(Fs/F0)/NSamples; % time to end SAG in [%]

n = 1:NSamples;
th_gmi_v(1:NSamples) = th_gmi; %threshold for magnitude detection
th_fi_v(1:NSamples) = th_fi; %threshold for frequency detection

Signal = SigGEN(F0,F1,Fs,Ps,NCycles,tau1,tau2,SNR,KaS, KxS);
%SigGEN(F0,F1,SampleRate,Ps,NCycles,tau1,tau2,SNR, KaS, KxS)

%calculate instantaneous freq and magnitude by Hilbert Transform
z=hilbert(Signal');  % calculates the analytic signal associated with Signal
fi = unwrap(angle(z)); gmi = abs(z);

% figure(1)
% plot(Signal,'.k'); ylabel('Sampled signal x[n] [V]'); xlabel('Samples')
% axis([1 480 -1 1]); grid on;

% subplot(3,1,1); plot(Signal); title('Signal'); grid on
% subplot(3,1,2); plot(fi); title('Instantaneous phase'); xlabel('Samples');ylabel('\theta_i[n]')
% subplot(3,1,3); plot(gmi); title('Instantaneous magnitude'); xlabel('Samples');ylabel('\a_i[n]')


%fazer o denoising por PATV
%% Perform PATV filtering
% PATV: Least square polynomial approximation + total variation denoising

% parameters
d_fi = 1;                          % d : degree of approximation polynomial
d_gmi = 0;
lam_fi = 1;
lam_gmi = 0.5;   % lambda : regularization parameter
Nit = 10;                      % Nit : number of iterations

%loPATV params
% a = 1.0;
% % a = eps;    % L1 norm - reproduces above result
% phi_fi = @(x) lam_fi/a * log(1 + a*abs(x));
% phi_gmi = @(x) lam_gmi/a * log(1 + a*abs(x));
% %dphi = @(x) lam ./(1 + a*abs(x)) .* sign(x);
% wfun_fi = @(x) abs(x) .* (1 + a*abs(x)) / lam_fi;
% wfun_gmi = @(x) abs(x) .* (1 + a*abs(x)) / lam_gmi;

%PATV algorithm
[x_fi, p_fi, cost_fi, u, v] = patv_MM(fi, d_fi, lam_fi, Nit);
[x_gmi,p_gmi, cost_gmi, u, v] = patv_MM(gmi, d_gmi, lam_gmi, Nit);

%loPATV algorithm
% [x_fi, p_fi, cost_fi, u, v] = patv_MM2(fi, d_fi,phi_fi,wfun_fi, Nit);
% [x_gmi,p_gmi, cost_gmi, u, v] = patv_MM2(gmi, d_gmi,phi_gmi,wfun_gmi, Nit);

% display cost function history
% figure(2)
% subplot(2,1,1)
% semilogy(1:Nit, cost_fi)
% title('PATV algorithm - Cost function history');
% xlabel('Iteration')
% ylabel('Cost function value')
% subplot(2,1,2)
% semilogy(1:Nit, cost_gmi)
% title('PATV algorithm - Cost function history');
% xlabel('Iteration')
% ylabel('Cost function value')

%display TV component
% figure(3)
% clf
% % subplot(3,1,1)
% % plot(n, fi,'.k', n, x_fi+p_fi, 'black') 
% % title('Calculated TV component (PATV) - Frequency');
% % legend('Data','Estimated signal', 'Location','southeast')
% subplot(2,1,1)
% plot(n, gmi,'.k', n, x_gmi+p_gmi, 'black') 
% title('Calculated TV component (PATV) - Magnitude');
% legend('Data','Estimated signal', 'Location','southeast')
detector_gmi = br_mask'.*abs(gradient(x_gmi));
% subplot(2,1,2)
% plot(n,detector_gmi,'b');  ylabel('Detection signal')
d2_gmi = abs(detector_gmi.*gradient(detector_gmi)); 

%hold on; plot(n,d2_gmi)

%%%%%%%%%%%  estimate taus

%detector_gmi
[gmi_max(1),imax_gmi(1)] = maxk(detector_gmi,1);
%[gmi_min,imin] = min(detector_gmi);
% ignores adjacent samples
%br = floor(0.02*NSamples);
aind = floor(samples_cycle)/2; %ignores half cycle centered in imax_gmi
detector_gmi(imax_gmi(1)-aind:imax_gmi(1)+aind) = 0;

% hold on; plot(n,detector_gmi,'red',n,th_gmi_v,'k--')
% legend('d1_{gmi}','d2_{gmi}', 'Location','southeast')
[gmi_max(2),imax_gmi(2)] = maxk(detector_gmi,1);

%detector_fi
detector_fi = br_mask'.*abs(gradient(x_fi));

% figure(4)
% subplot(2,1,1)
% plot(n, fi,'.k', n, x_fi+p_fi, 'black') 
% title('Calculated TV component (PATV) - Frequency');
% legend('Data','Estimated signal', 'Location','southeast')
% subplot(2,1,2)
% plot(n,detector_fi,'b');  ylabel('Detection signal')

[fi_max(1),imax_fi(1)] = maxk(detector_fi,1); %detects first peak

%ignore adjacent samples
%detector_fi(imax_fi-br:imax_fi+br) = 0; 
aind = imax_fi-floor((Fs/F0)/NSamples/2);
while detector_fi(aind)>th_fi
    detector_fi(aind:imax_fi) = 0;
    aind = aind - 1;
end
aind = imax_fi+1;
while detector_fi(aind)>th_fi
    detector_fi(imax_fi:aind) = 0;
    aind = aind + 1;
end
% hold on; plot(n,detector_fi,'red',n,th_fi_v,'k--')
% legend('d1_{fi}','d2_{fi}', 'Location','southeast')
[fi_max(2),imax_fi(2)] = maxk(detector_fi,1); %detects second peak

% error calculation
[tau_est_gmi, gmi_ind] = sort(imax_gmi); %error in samples
tau_est_gmi(gmi_max < th_gmi) = NaN; % check if estimates are valid
tau_error_gmi = (tau_est_gmi - [tau1 tau2]*NSamples); %error in [dt]

[tau_est_fi,fi_ind] = sort(imax_fi); %tau estimated in [samples]
tau_est_fi(fi_max < th_fi) = NaN; % check if estimates are valid
tau_error_fi = (tau_est_fi - [tau1 tau2]*NSamples); %error in [dt]

%criterio de escolha mag x freq com denoising: usar a maior diferença para 
% o limiar
crit = gmi_max(gmi_ind) - th_gmi > fi_max(fi_ind) - th_fi;
tau_est(crit) = tau_est_gmi(crit);
tau_est(~crit) = tau_est_fi(~crit);
tau_error = tau_est - [tau1 tau2]*NSamples;

%calcular frequencia pela inclinação da fase
P = polyfit(n,p_fi',1);
f_est = P(1)*Fs/(2*pi);
FE = (f_est - F1)*100/F1;

tau_est = tau_est/Fs;

%problemas: 
% 1 - indices excedendo as bordas (na hora de excluir do sinal de
% detecção) -> implementar um teste para isso não ocorrer
% por exemplo, baseado no th_gmi - ok, mas usando meio ciclo em torno do
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
