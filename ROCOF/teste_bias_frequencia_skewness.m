%testando o bias na estimacao da frequencia
clear all; close all; clc;

%fixed parameters
SNR = 90;
F0 = 60.0; %nominal frequency
F1 = 60.; %fundamental frequency
KaS = 0.0; %[degrees]
KxS = 0.0; % [relative magnitude step]
KfS = 0.0; %[Hz] %size of frequency step
Fs =4800;
NCycles = 6;
T = NCycles/F0;
phi_0 = 70; %angle phi_0 [degrees]
tau1 = 0.5;  % in proportion of T
tau2 = 1.; % in proportion of T; set tau2 = 1 if you dont want two steps
NSamples = floor(NCycles*Fs/F0);
tau_n1 = floor(tau1*NSamples);
tau_n2 = floor(tau2*NSamples);
nbits = 16;
br = 80/480; brn = floor(br*NSamples); % number of samples to be ignored 
brmask = [(brn+1):(tau_n1-brn+1) ((tau_n1+brn+2):(NSamples-brn))];
%brmask = [(brn+1:tau_n1-brn+1) (tau_n1+brn:NSamples-brn+2)];

for k = 1:5
Signal = SigGEN2(F0,F1,Fs,phi_0,NCycles,tau1,tau2,SNR,KaS, KxS,KfS,nbits);
Sfft(k,:) = fft(Signal');
z=hilbert(Signal');

Psi_i = unwrap(angle(z)); %[rad]
%Psi_i = phase(z);
f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
az = abs(z);
f_i1 = f_i.*az./median(az);
sk1(k) = skewness(f_i1);
SNR = SNR - 10;
%KfS = KfS + 1;
%KaS = KaS + 10;
fest1 = median(f_i1(brn:(end-brn)));
FE1(k) = fest1 - F1;

%brmask = [(tau_n1+brn+2:NSamples-brn)]; %estima a freq somente após o tau1
az2 = az(brmask);
f_i2 = f_i(brmask);
%f_i2 = [f_i(1:tau_n1-brn+1); f_i(tau_n1+brn:end)];
fest2 = median(f_i2)
f_i3 = f_i2.*az2./median(az2);
%plot(f_i2,'.'); hold on;
fest3 = median(f_i3)
FE2(k) = fest2 - F1;
FE3(k) = fest3 - F1;
end

%legend('KaS = -10', 'KaS = 0', 'KaS = +10')

%figure; plot(f_i2); hold on; plot(f_i3); legend('f_i2','f_i3')

FE1
FE2
FE3

%sk = skewness(f_i)
%sk1 = skewness(f_i1)
%sk2 = skewness(f_i2)
%sk3 = skewness(f_i3)

% figure
% df = Fs/NSamples;
% fx = 0:df:Fs/2-df;
% semilogy(fx,abs(Sfft(1,1:NSamples/2))); 
% hold on; semilogy(fx,abs(Sfft(3,1:NSamples/2)))
% legend('FFt KaS = -10','FFT KaS = 0')