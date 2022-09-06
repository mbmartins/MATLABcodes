%generation of performance comparison
% avaliação da distribuição de ri
clear all; close all; clc

SNR = 40:10:100;

%angulo phi_0, fixo ou tabela
Ps = 0;  %pior caso para salto de freq é phi0 = 0;

%fixed parameters
F0 = 60.0;
F1 = 60.0;
Fs = 4800;
NCycles = 6;
tau = 0.5;  % in [%] of the time window
tau2 = 10;
SAG_cycles = 10; %duration of SAG 
nbits = 16;

% %CASE 3 - salto de frequencia
h_a = 0.0; %[degrees]
h_x = -0.; % [relative step]
h_f = -1.; %[Hz]

for s = 1:length(SNR)  % loop for different SNRs
    rng default; 
    Signal = SigGEN2(F0,F1,Fs,Ps,NCycles,tau,tau2,SNR(s),h_a,h_x,h_f,nbits);

    N = length(Signal);
    br = floor(0.1*N); % 5% of NSamples are taken off at the beggining and end
    %excluindo amostras do inicio e do final
    brmask = br:N-br-1;

    z=hilbert(Signal');  % calculates the analytic signal associated with Signal
    wi=gradient(unwrap(angle(z)));% Hilbert estimate of the instantaneous frequency of z
    ai = abs(z);
    fi(s,:) = wi(brmask);

    ri(s,:) = gradient(fi(s,:)); %rocof de fi compensado
    %pd = fitdist(ri,'Normal');
    [h(s),p(s)] = chi2gof(ri(s,:));
end
n = 1:N;
nn = n(brmask);
% figure(1);
% plot(nn,ri(3,:)); hold on;
% plot(nn,ri(5,:));
% legend('SNR = 60 dB','SNR = 80 dB')
tau_n = N*tau;
degrau = [ones(1,N)]; degrau(tau_n:N) = (F0+h_f)/F0;
f_degrau = 60*degrau;
f_degrau = f_degrau(brmask);
r_degrau = gradient(f_degrau);

figure; 
plot(nn,f_degrau,'k--'); hold on;
plot(nn,Fs*fi(end,:)/(2*pi)); 
plot(nn,Fs*fi(3,:)/(2*pi),'r:'); 
ylabel('$f_i[n]$ [Hz]','Interpreter','latex','Fontsize',14)
xlabel('Amostras')
lg= legend('$f_i[n]$','$\hat{f}_i[n]$, SNR = 100 dB', '$\hat{f}_i[n]$, SNR = 80 dB')
lg.Interpreter = 'latex';

figure; 
plot(nn,abs(r_degrau),'k--'); hold on;
plot(nn,abs(Fs*ri(end,:)/(2*pi)),'LineWidth',2); 
plot(nn,abs(Fs*ri(5,:)/(2*pi)),'r:','LineWidth',1); 
ylabel('$|r_i[n]|$ [Hz/s]','Interpreter','latex','Fontsize',14)
xlabel('Amostras')
lg2= legend('$r_i[n]$','$\hat{r}_i[n]$, SNR = 100 dB', '$\hat{r}_i[n]$, SNR = 80 dB')
lg2.Interpreter = 'latex';
lg2.FontSize = 13;
lg2.Location = 'Northwest';
grid on;