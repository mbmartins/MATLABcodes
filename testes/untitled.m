clear all; close all; clc;

f0 = 50;
w0 = 2*pi*f0;
fm = 5;
w = 2*pi*fm;
dt = 1/5e3;
t = 0:dt:0.1;
phi=0;

ka = 1; %há um limite máximo para ka?
dka = -5:0.1:5;
dw = -2*pi:0.1:pi;

for j = 1:length(dw);
for i = 1:length(dka);

ka_est = ka + dka(i);
w_est = w + dw(j);

f = cos(w0*t + ka*sin(w*t + phi));
f_est = cos(w0*t + ka_est*sin(w_est*t + phi));

Eps(i,j) = sum((f - f_est).^2);
end
end

surf(Eps)
xlabel('ka'); ylabel('w'); zlabel('E^2')