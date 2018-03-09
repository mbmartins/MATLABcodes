% Estimation of Sincrophasors with discontinuities
clear all; close all; clc;

%signal generation
F0 = 60; %Hz
F1 = 60; %Hz
SampleRate = 4800; %Hz
dt = 1/SampleRate;
AnalysisCycles = 4;
NSamples = floor(AnalysisCycles*SampleRate/F0);
n = -NSamples/2:(NSamples/2-1); %discrete time vector
tau_0 = 0; %discrete time displacement
n = n + tau_0;
t = n*dt; %time vector
Vm = 100; %70*sqrt(2);
Ps = 0; %phase in degrees

% Phase in radians
Ph = Ps*pi/180;

KaS = 30;   % phase (angle) step index
KxS = 0.0;   % magnitude step index
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
% Phase step
% Model: y = B*cos(w*t + phi + A*u(t - tau))

%first estimates for A, B, and tau
A = 30;
Vm = 100;
B = Vm;
tau = 0;

%estimated signal
u = zeros(length(Xm),length(t));
u(length(Xm),t >= tau) = u(length(Xm),t >= tau) + 1;
y = B*cos(Wf*t + Ph + A*(pi/180)*u);

plot(t,Signal,'.b',t,y,'r')

%sum of squared errors with variation of A and B

%using A as parameter
dA = 0.1;
A = (27:dA:33);

dB = 0.001;
B = (0.980:dB:1.01)*Vm;
%dtau = dt;
%tau = -50*dtau:dtau:49*dtau;

for i = 1:length(A)   
%    for j = 1:length(tau)
for j = 1:length(B)
        %u = zeros(length(Xm),length(t));
        %u(length(Xm),t >= tau(j)) = u(length(Xm),t >= tau(j)) + 1;
        y = B(j).*cos(Wf*t + Ph + A(i)*(pi/180)*u);
        err(i,j) = sum((Signal - y).^2);
end
%    end
end

figure
%plot(A,err)
surf(B,A,err)
title('Error surface: A x B x error')
% 
% %gradient descendent
% dedA = sum((y - B*(1+A*u)*cos(
% dB
% dtau
B = 101;
A = 31;
y = B.*cos(Wf*t + Ph + A*(pi/180)*u);
figure
plot(t,Signal - y,'r')
title('Point by point estimation error: Signal - y')

eps = 0.001;
err = sum((Signal - y).^2);
while (err)>eps
    y_ = B.*cos(Wf*t + Ph + A*(pi/180)*u);
    dedB = sum((Signal - y_).*cos(Wf*t + Ph + A*(pi/180)*u));
    dedA = sum(-(Signal - y_).*sin(Wf*t + Ph + A*(pi/180)*u)*(pi/180).*u);
    
    %finding the best alfa:
    alfa1 = 0.1;
    alfa2 = 0.01;
%     p = 0.01;
%     derrdalfa = sum(Signal - B*cos(Wf*t + Ph + (A+alfa*dedA)*(pi/180)*u))
%     alfa = 
    A = A + alfa1*dedA
    B = B + alfa2*dedB
    err = sum((Signal - y_).^2)
end    
    