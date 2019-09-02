%% sim PLL

f1 = 60.;
ef = 0.00;
phi = 60;

Fs = 5000;
dt = 1/Fs;
t = (1:5000)*dt;

ruido = 1.e-3*randn(1,length(t));
signal = cos(2*pi*(f1+ef)*t + phi*pi/180) + ruido;

%mede freq por hilbert
z = hilbert(signal);
phi_i = unwrap(angle(z));
%plot(t,phi_i)
fi = gradient(phi_i)/dt/(2*pi);
fmed = mean(fi)

%erro_fmedido = 0.e-4;
%fmedido = f1+ef + erro_fmedido;
fmedido = f1

x = signal*cos(2*pi*fmedido*t)';
y = signal*sin(2*pi*fmedido*t)';

fasor = x - j*y;

fase = angle(fasor)*180/pi

FE = fmedido - f1
erro_fase = fase - phi

