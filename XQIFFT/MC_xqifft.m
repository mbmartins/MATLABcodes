clear all; close all; clc;

%TODO :ver influencia em p e FE da variação em f1
% verificar a geração do ruido se esta correta
% comparar a medição com degraus

%signal generation
N = 2.^(8:11); %number of samples
%cycles = 1;
fs = 2500; %sampling frequency;
f0 = 50.; %nominal system frequency
Uf = 0.01;

for ni = 1:size(N,2)
    dt = 1/fs;
    fbin = fs/N(ni);
%Monte Carlo 
for h = 1:1
    f1 = f0 + Uf*(randn(1,1)-0.5); %actual fund frequency
    n = (0:N(ni)-1);
    t = n*dt;
    u = [zeros(1,N(ni)/2) ones(1,N(ni)/2)];
    kx = 0.0; ka = 0*pi/180;%[rad]
    Vm = 1;
    SNR = 60; %[dB]
    Aruido = Vm/10^(SNR/20);
    x = Vm*(1+kx*u).*sin(2*pi*f1*t+pi/2+ka*u); %+ 0.5*sin(2*pi*36*f1*t);  %samples
    x = x +  Aruido*(randn(1,length(t))-0.5); % ???????
    %plot(t,x,'o--')
    %snr_sig = snr(x)

    %windowing
    %rectangular
        window_rect = ones(1,N(ni));
    %hanning
        window_hanning = (0.5 - 0.5*cos(2*pi*n/N(ni)));
        % pmin = 0.246;
    %blackman-harris
        window_bm = blackmanharris(N(ni))';
        %wvtool(window);
        % pmin = 0.0853

    ps = 1e-4; pini = 0.05; pfin = 1;
    p = pini:ps:pfin;
    %p = 0.2240;
    for k = 1:size(p,2)
        [f_ip(k),A_ip(k),ph_ip(k)] = ipfft(x,fs,p(k),window_bm,'parabola');
        [f_lq(k),A_lq(k),ph_lq(k)] = ipfft(x,fs,p(k),window_bm,'log');
        [f_xq_bm(k),A_xq_bm(k),ph_xq_bm(k)] = ipfft(x,fs,p(k),window_bm,'power');
        [f_xq_han(k),A_xq_han(k),ph_xq_han(k)] = ipfft(x,fs,p(k),window_hanning,'power');
    end

    fe_ip = (f_ip - f1);
    fe_lq = (f_lq - f1);
    fe_xq_bm = (f_xq_bm - f1);
    fe_xq_han = (f_xq_han - f1);

    [fe_xqbm_min(h),ipmin] = min(abs(fe_xq_bm));
    pmin_bm(h,ni) = p(ipmin);
    fest_bm = f_xq_bm(ipmin);
    fe_bm(h,ni) = fest_bm - f1;

    [fe_xqhan_min(h),ipmin] = min(abs(fe_xq_han));
    pmin_han(h,ni) = p(ipmin);
    fest_han = f_xq_han(ipmin);
    fe_han(h,ni) = fest_han - f1;
    h
end
end
% figure
% subplot(2,1,1)
% hist(fe_han); xlabel('FE_{HAN}')
% mean(fe_han)
% std(fe_han)
% subplot(2,1,2)
% hist(pmin_han); xlabel('p value')
% mean(pmin_han)
% std(pmin_han)

% figure
% subplot(2,1,1)
% hist(fe_bm); xlabel('FE_{BM}')
% mean(fe_bm)
% std(fe_bm)
% subplot(2,1,2)
% hist(pmin_bm); xlabel('p value')
% mean(pmin_bm)
% std(pmin_bm)


figure
loglog(p,abs(fe_ip),p,abs(fe_lq),p,abs(fe_xq_bm),p,abs(fe_xq_han)); 
xlabel('p value'); ylabel('fe [Hz]')
xlim([pini pfin]);
legend('IP-DFT','LQ-DFT','XQ-DFT-BM','XQ-DFT-HAN')

figure
plot(N,pmin_bm,'o--',N,pmin_han,'x--')
xlabel('N amostras'); ylabel('p value')
legend('XQ-DFT-BM','XQ-DFT-HAN')

figure
semilogy(N,abs(fe_bm),'o--',N,abs(fe_han),'x--')
xlabel('N amostras'); ylabel('FE[Hz]')
legend('XQ-DFT-BM','XQ-DFT-HAN')

% o valor p depende principalmente de:
% 1 - qual a janela utilizada: hanning, blackmann-harris, etc
% 2 - numero de amostras

% numero de ciclos, 
% tamanho do afundamento kx
% tamanho do degrau de fase ka
% deltaf = f0 - fest