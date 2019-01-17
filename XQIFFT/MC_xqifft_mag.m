clear all; close all; clc;

%TODO :ver influencia em p e FE da variação em f1
% verificar a geração do ruido se esta correta
% comparar a medição com degraus

%signal generation
N = 512; %number of samples
%cycles = 1;
fs = 5000; %sampling frequency;
f0 = 50.0; %nominal system frequency
Uf = 0.0;
Vm = 1.;
    ni = 1;
dt = 1/fs;
    fbin = fs/N(ni);

    
%Monte Carlo 
for h = 1:1
    f1 = f0 + Uf*(randn(1,1)-0.5); %actual fund frequency
    n = (0:N(ni)-1);
    t = n*dt;
    u = [zeros(1,N(ni)/2) ones(1,N(ni)/2)];
    kx = 0.0; ka = 0*pi/180;%[rad]
    x = Vm*(1+kx*u).*sin(2*pi*f1*t+ka*u); %+ 0.5*sin(2*pi*36*f1*t);  %samples
    var_sig = std(x);
    SNR = 60; %[dB]
    eta = var_sig/10^(SNR/20); %eq (3) CPEM
    x = x +  eta*(randn(1,length(t))-0.5); 
    %plot(t,x,'o--')
    snr_sig = snr(x)

    %windowing
    %rectangular
        window_rect = ones(1,N(ni));
    %hanning
        window_hanning = (0.5 - 0.5*cos(2*pi*n/N(ni)));
        % pmin = 0.246; para freq
    %blackman-harris
        window_bm = blackmanharris(N(ni))';
        %wvtool(window);
        % pmin = 0.0853 para freq

    ps = 1e-4; pini = 0.05; pfin = 1;
    p = pini:ps:pfin;
    %p = 0.2240;
    for k = 1:size(p,2)
        [f_ip(k),A_ip(k),ph_ip(k)] = ipfft(x,fs,p(k),window_bm,'parabola');
        [f_lq(k),A_lq(k),ph_lq(k)] = ipfft(x,fs,p(k),window_bm,'log');
        [f_xq_bm(k),A_xq_bm(k),ph_xq_bm(k)] = ipfft(x,fs,p(k),window_bm,'power');
        [f_xq_han(k),A_xq_han(k),ph_xq_han(k)] = ipfft(x,fs,p(k),window_hanning,'power');
    end

    Ae_ip = (A_ip - Vm);
    Ae_lq = (A_lq - Vm);
    Ae_xq_bm = (A_xq_bm - Vm);
    Ae_xq_han = (A_xq_han - Vm);

    [Ae_xqbm_min(h),ipmin] = min(abs(Ae_xq_bm));
    pmin_bm(h,ni) = p(ipmin);
    Aest_bm = A_xq_bm(ipmin);
    Ae_bm(h,ni) = Aest_bm - Vm;

    [Ae_xqhan_min(h),ipmin] = min(abs(Ae_xq_han));
    pmin_han(h,ni) = p(ipmin);
    Aest_han = A_xq_han(ipmin);
    Ae_han(h,ni) = Aest_han - Vm;
    h
end

% figure
% title('AE IPDFT')
% hist(Ae_ip)
% 
% figure
% title('AE LQ-DFT')
% hist(Ae_lq)

figure
subplot(2,1,1)
hist(Ae_han); xlabel('AE_{HAN}')
mean(Ae_han)
std(Ae_han)
subplot(2,1,2)
hist(pmin_han); xlabel('p value')
mean(pmin_han)
std(pmin_han)

figure
subplot(2,1,1)
hist(Ae_bm); xlabel('AE_{BM}')
mean(Ae_bm)
std(Ae_bm)
subplot(2,1,2)
hist(pmin_bm); xlabel('p value')
mean(pmin_bm)
std(pmin_bm)

figure
loglog(p,abs(Ae_ip),p,abs(Ae_lq),p,abs(Ae_xq_bm),p,abs(Ae_xq_han)); 
xlabel('p value'); ylabel('|Ae| [V]')
xlim([pini pfin]);
legend('IP-DFT','LQ-DFT','XQ-DFT-BM','XQ-DFT-HAN')

% o valor p depende principalmente de:
% 1 - qual a janela utilizada: hanning, blackmann-harris, etc
% 2 - numero de amostras

%investigação sobre a fase seria nova!

% numero de ciclos, 
% tamanho do afundamento kx
% tamanho do degrau de fase ka
% deltaf = f0 - fest