clear all; close all; clc;

%testando a aplicação de PATV para o sinal de magnitude instantânea

load('IM_sig.mat')
y = IM_sig;
n = 1:length(y);

figure(1)
plot(n, y, '.k')
title('Noisy data');

%% Perform PATV filtering
% PATV: Least square polynomial approximation + total variation denoising

% parameters
d = 1;                          % d : degree of approximation polynomial
lam = 1;                      % lambda : regularization parameter
Nit = 100;                      % Nit : number of iterations

[x, p, cost, u, v] = patv_MM(y, d, lam, Nit);

% display cost function history
figure(2)
plot(1:Nit, cost)
title('PATV algorithm - Cost function history');
xlabel('Iteration')
ylabel('Cost function value')
%ylim([3.3 4.5])
%printme('convergence')

%%
% Display calculated TV component

figure(3)
clf
subplot(2,1,1)
plot(n, x, 'k')
title('Calculated TV component (PATV)');
%printme('TV_component')

detector = gradient(x);
subplot(2,1,2)
plot(n,detector)

[dmax,imax] = max(detector)



%%
% Display TV-compensated data

figure(4)
clf
plot(n, y - x, '.k')
title('TV-compensated data (PATV)');
%printme('TV_compensated_data')


%%
% Display estimated signal

figure(5)
clf
plot(n, y,'.k', n, x+p, 'black') 
title('Estimated signal (PATV)');
legend('Data','Estimated signal', 'Location','south')
%printme('PATV')