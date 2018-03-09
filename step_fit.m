% Estimation of Sincrophasors with discontinuities
clear all; close all; clc;

%signal generation
F0 = 60; %Hz
F1 = 60; %Hz
SampleRate = 4800; %Hz
dt = 1/SampleRate;
AnalysisCycles = 40;
NSamples = floor(AnalysisCycles*SampleRate/F0);
n = -NSamples/2:(NSamples/2-1); %discrete time vector
t = n*dt; %time vector
Vm = 100; %70*sqrt(2);
Ps = 0; %phase in degrees

% Phase in radians
Ph = Ps*pi/180;

KaS = 0.0;   % phase (angle) step index
KxS = 0.1;   % magnitude step index
Wf = 2*pi*F1;  % fundamental frequency

Xm = Vm; %for now, single phase; TODO: 6-channels
Ain = zeros(length(Xm),length(t));
% Amplitude Step: applied after time passes 0
i = 1;
Ain(i,:) = Xm(i);
Ain(i,t >= 0) = Ain(i,t >= 0) * (1 + KxS(i));
%Phase step
Theta(i,:) = (Wf(i)*t) ...                         % Fundamental
                 + Ph(i);               % phase shift
Theta(i,t >= 0) = Theta(i,t >= 0) + (KaS(i) * pi/180);
cSignal = (Ain.*exp(-1i.*Theta));
Signal = real(cSignal);

%%%% Signal parameters estimation
% Amplitude step
% Model: y = B*(1 + A*u(t - tau)*cos(w*t + phi))

%first estimates for A, B, and tau
A = KxS*0.99;
Vm = 100;
B = Vm;
tau = 0;

%estimated signal
u = zeros(length(Xm),length(t));
u(length(Xm),t >= tau) = u(length(Xm),t >= tau) + 1;
y_ = B*(1 + A*u).*cos(Wf*t + Ph);

plot(t,Signal,'b',t,y_,'r')

%sum of squared errors with variation of A and tau 

%using A/B as parameter
dA = 0.0005;
A = (0.09:dA:0.11);

dB = 0.001;
B = (0.990:dB:1.01)*Vm;
%dtau = dt;
%tau = -50*dtau:dtau:49*dtau;

for i = 1:length(A)   
%    for j = 1:length(tau)
for j = 1:length(B)
        %u = zeros(length(Xm),length(t));
        %u(length(Xm),t >= tau(j)) = u(length(Xm),t >= tau(j)) + 1;
        y_ = B(j)*(1 + A(i)*u).*cos(Wf*t + Ph);
        dedA(i,j) = sum((Signal - B(j)*(1+A(i)*u).*cos(Wf*t + Ph)).*(-B(j)*u.*cos(Wf*t + Ph)));
        err(i,j) = sum((Signal - y_).^2);
end
%    end
end

figure
%plot(A,err)
% subplot(1,2,1)
surf(B,A,err);
title('Error surface: A x B x error')
% subplot(1,2,2)
% surf(B,A,dedA);
% title('Error derivative: dedA x error')

% 
% %Nonlinear gradient descendent
% A = 0.11; B = 99;
% dedA = sum((y - B*(1+A*u).*cos(Wf*t + Ph)).*(B*u.*cos(Wf*t + Ph)))
% dedB = sum((y - B*(1+A*u).*cos(Wf*t + Ph)).*(1 + A*u.*cos(Wf*t + Ph)))
% dtau

eps = 1; %stop condition
A = 0.1; 
B = 99.9;
alfa = 1e-7;
beta = 1e-3;
err = 1000;

while (err)>eps
    y_ = B*(1 + A*u).*cos(Wf*t + Ph);
    dedA = sum((Signal - B*(1+A*u).*cos(Wf*t + Ph)).*(-B*u.*cos(Wf*t + Ph)));
    dedB = sum((Signal - B*(1+A*u).*cos(Wf*t + Ph)).*(-1-A*u.*cos(Wf*t + Ph)));
    dx = -[dedA dedB];
    %line search alfa
    %alfa = argmin e(x + alfa*dx)
    dedalfa = 1; tol = 1e-3;
    while dedalfa>tol
        alfa = sqrt(dx(1)^2+dx(2)^2)/1e6; %initial guess
        A_ = A + alfa*dx(1);
        B_ = B + alfa*dx(2);
        e = sum((Signal - B_*(1 + A_*u).*cos(Wf*t + Ph)).^2);
        dalfa = alfa/100;
        A_ = A_ + (alfa + dalfa)*dx(1);
        B_ = B_ + (alfa + dalfa)*dx(2);
        e_ = sum((Signal - B_*(1 + A_*u).*cos(Wf*t + Ph)).^2);
        dedalfa = (e_ - e)/dalfa; %gradient descent for alfa
    end
    
    [dA dB]= + alfa*dx
    [A B] = [A + dA B + dB]
    err = sum((Signal - y_).^2)
end


% figure
% plot(t,Signal - y_,'r')
% title('Point by point estimation error: Signal - y')


