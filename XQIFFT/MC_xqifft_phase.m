clear all; close all; clc;

% Verificação da implementação de XQIFFT
% Variáveis de influência na estimação de Frequência
% Estudo da influência da variação de p em FE 
% usando p fixo para cada realização de MC
% deltaf = f0 - fest

%signal generation
%cycles = 1;
fs = 5000; %sampling frequency;
N = 5000; %number of samples
f0 = 50.0; %nominal system frequency
Uf = 0.005; %uncertainty of frequency
Vm = 1.;
    ni = 1;
dt = 1/fs;
    fbin = fs/N(ni);
    
f1 = f0 + (-50:50)*Uf*f0/100;
%deltaf = f1 - f0;

for ni = 1:size(N,2)
    dt = 1/fs;
    fbin = fs/N(ni);
%Monte Carlo 
for h = 1:101
    %f1(h) = f0 + Uf*f0*(randn(1,1)-0.5); %actual fund frequency
    n = (0:N(ni)-1);
    t = n*dt;
    u = [zeros(1,N(ni)/2) ones(1,N(ni)/2)];
    kx = 0.0; ka = 0*pi/180;%[rad]
    x = Vm*(1+kx*u).*sin(2*pi*f1(h)*t+ka*u + pi/2); %+ 0.5*sin(2*pi*36*f1*t);  %samples
    var_sig = std(x);
    SNR = 60; %[dB]
    eta = var_sig/10^(SNR/20); %eq (3) CPEM
    x = x +  eta*(randn(1,length(t))-0.5);
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

    %Busca do valor p que minimiza FE
%     ps = 1e-4; pini = 0.05; pfin = 1;
%     p = pini:ps:pfin;
    
    p = 0.2240;
    for k = 1:size(p,2)
        [f_ip(k),A_ip(k),ph_ip(k,h)] = ipfft(x,fs,p(k),window_hanning,'parabola');
        [f_lq(k),A_lq(k),ph_lq(k,h)] = ipfft(x,fs,p(k),window_hanning,'log');
        [f_xq_bm(k),A_xq_bm(k),ph_xq_bm(k,h)] = ipfft(x,fs,p(k),window_bm,'power');
        [f_xq_han(k),A_xq_han(k),ph_xq_han(k,h)] = ipfft(x,fs,p(k),window_hanning,'power');
    end

    fe_ip(h) = (f_ip - f1(h));
    fe_lq(h) = (f_lq - f1(h));
    fe_xq_bm = (f_xq_bm - f1(h));
    fe_xq_han = (f_xq_han - f1(h));

    [fe_xqbm_min(h),ipmin] = min(abs(fe_xq_bm));
    pmin_bm(h,ni) = p(ipmin);
    fest_bm = f_xq_bm(ipmin);
    fe_bm(h,ni) = fest_bm - f1(h);

    [fe_xqhan_min(h),ipmin] = min(abs(fe_xq_han));
    pmin_han(h,ni) = p(ipmin);
    fest_han(h) = f_xq_han(ipmin);
    fe_han(h,ni) = fest_han(h) - f1(h);
    h
end
end

%figure; plot(f1,ph_xq_han*180/pi,'o');xlabel('F_1 [Hz]');ylabel('phase [deg]');
%figure; plot(f1,ph_ip*180/pi,'o');xlabel('F_1');ylabel('phase [deg]');

figure; plot(f1,ph_lq,'o');xlabel('F_1');ylabel('phase [rad]');
ph_corr = ph_lq - pi*(fest_han-f0);
figure;plot(f1,ph_corr)

% figure
% % subplot(2,1,1)
% hist(fe_han); xlabel('FE_{HAN}')
% title('FE_{HAN}, p = 0.2240')
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

% 
% figure
% loglog(p,abs(fe_ip),p,abs(fe_lq),p,abs(fe_xq_bm),p,abs(fe_xq_han)); 
% xlabel('p value'); ylabel('FE [Hz]')
% xlim([pini pfin]);
% legend('IP-DFT','LQ-DFT','XQ-DFT-BM','XQ-DFT-HAN')

% figure
% plot(N,pmin_bm,'o--',N,pmin_han,'x--')
% xlabel('N amostras'); ylabel('p value')
% legend('XQ-DFT-BM','XQ-DFT-HAN')
% 
% figure
% semilogy(N,abs(fe_bm),'o--',N,abs(fe_han),'x--')
% xlabel('N amostras'); ylabel('FE[Hz]')
% legend('XQ-DFT-BM','XQ-DFT-HAN')

% COMENTARIOS
% o valor p depende principalmente de:
% 1 - qual a janela utilizada: hanning, blackmann-harris, etc
% 2 - numero de amostras

% numero de ciclos, 
% tamanho do afundamento kx
% tamanho do degrau de fase ka
% deltaf = f0 - fest