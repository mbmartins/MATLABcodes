%teste PATV
clc; close all; clear all;

a = 1.5; b = 2;
t = [1:10];
var = 0.1;
r = var*randn(1,10);
deg = [0 0 0 0 0 1 1 1 1 1];
pol = a*t + b;
polr = pol + r + deg;
plot(polr)
hold on;
plot(pol); plot(deg);

d = 1; lambda = 1.e-5; Nit = 10;
[x, p, cost, u, v] = patv_MM(polr, d, lambda, Nit);

figure; hold on;
plot(x)
%plot(p)
plot(u)
%plot(v)
legend('x','u','v')

pol_est = p + x;

figure; hold on;
plot(t,polr,'b',t,pol_est,'r')
legend('pol', 'p estimado')

deg_max = max(deg);
u_max = max(u);

deg_error = deg_max - u_max