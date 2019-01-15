clear all; close all; clc;

%Influencia de delta_fbin, usando hanning window

%signal generation
fnom = 60.0; %nominal fund frequency
fs = 5000; %sampling frequency;
dt = 1/fs;
N = 500; %number of samples
fbin = fs/N;

for fk=1:1000
    f1(fk) = fnom + 2*(fk/1e3)*fbin - fbin;
    n = (0:N-1);
    t = n*dt;
    u = [zeros(1,N/2) ones(1,N/2)];
    kx = 0.0; ka = 0*pi/180;%[rad]
    SNR = 60; %[dB]
    x = (1+kx*u).*sin(2*pi*f1(fk)*t+pi/2+ka*u) + 0.1*sin(2*pi*36*f1(fk)*t);  %samples
    x = x + awgn(x,SNR,'measured');
    %plot(t,x)

    %windowing
    %rectangular
%        window_rect = ones(1,N);
    %hanning
        window_hanning = (0.5 - 0.5*cos(2*pi*n/N));
        % pmin = 0.246;
    %blackman-harris
        window_bm = blackmanharris(N)';
        %wvtool(window);
        % pmin = 0.0853

    ps = 1e-4; pini = 0.05; pfin = 1;
    %p = pini:ps:pfin;
    phan = 0.2240;
    pbm = 0.0853
    [f_ip,A_ip,ph_ip] = ipfft(x,fs,phan,window_bm,'parabola');
    [f_lq,A_lq,ph_lq] = ipfft(x,fs,phan,window_bm,'log');
	[f_xq_bm,A_xq_bm,ph_xq_bm] = ipfft(x,fs,pbm,window_bm,'power');
	[f_xq_han,A_xq_han,ph_xq_han] = ipfft(x,fs,phan,window_hanning,'power');

    fe_ip = (f_ip - f1);
    fe_lq = (f_lq - f1);
    fe_xq_bm(fk) = (f_xq_bm - f1(fk));
    fe_xq_han(fk) = (f_xq_han - f1(fk));

end

plot(f1,fe_xq_bm,'r',f1,fe_xq_han,'b')
xlabel('Frequency [Hz]'); ylabel('FE [Hz]')
xlim([fnom-0.1 fnom+0.1])
legend('FE_{BM}','FE_{HAN}')

% figure
% loglog(p,fe_ip,p,fe_lq,p,fe_xq_bm,p,fe_xq_han); 
% xlabel('p value'); ylabel('fe [Hz]')
% xlim([pini pfin]);
% legend('IP-DFT','LQ-DFT','XQ-DFT-BM','XQ-DFT-HAN')

% o valor p depende principalmente de:
% qual a janela utilizada: hannin, blackmann-harris, etc
% numero de ciclos
% tamanho do afundamento kx
% tamanho do degrau de fase ka
% deltaf = f0 - fest

% estudo de p x deltaf
