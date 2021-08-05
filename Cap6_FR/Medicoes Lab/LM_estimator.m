function [params] = LM_estimator(Signal,x0,tau_n, KaS, KxS, Fs)
% x0 é o vetor com os valores iniciais dos parâmetros
% para salto de fase: xnom = [Vm 2*pi*F1 Ph KaS];
% para salto de magnitude: xnom = [Vm KxS 2*pi*F1 Ph];
x0(1) = max(Signal);
%tau_n é o índice no qual ocorre o salto em Signal
N = length(Signal);
tau_pp = tau_n/N; % relative time of step in percent of total time 
tau_0 = (tau_pp - 0.5)*N; %discrete time displacement
n = -N/2:(N/2-1); %discrete time vector
n = n - floor(tau_0);
dt = 1/Fs;
t = n*dt; %time vector

u = zeros(1,N);
tau = 0; %depois do deslocamento do vetor n, tau é sempre em zero
u(1,t >= tau) = u(1,t >= tau) + 1;

%modelo 1
% Phase step
% Model: f(x) = x1*cos(w*t + phi + x2*u(t - tau)) 
if KaS ~= 0
    %f = @(x) x(1)*cos(x(2)*t + x(3) + x(4)*(pi/180)*u);
    f = @(x) x(1)*cos(x(2)*t + x(3) + 2*pi + x(4)*(pi/180)*u);
else
    %mag
    %f = @(x) x(1)*(1+x(4)*u).*cos(x(2)*t + x(3));
    %xnom1 = [Vm 2*pi*F1 Ph KxS];
    f = @(x) x(1)*(1+x(2)*u).*cos(x(3)*t + x(4));
end
err = @(x) (Signal - f(x));

% Levenberg-marquardt from optimization toolbox
tol = 1e-7;
OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','OptimalityTolerance',tol);
OPTIONS.StepTolerance = 1e-12;
[X,RESNORM,RESIDUAL,exitflag,output] = lsqnonlin(err,x0,[],[],OPTIONS);
Y = f(X);

params = X;