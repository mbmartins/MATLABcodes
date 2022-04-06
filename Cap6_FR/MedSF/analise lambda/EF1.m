function [f1,f2,f_r,f_u,ri] = EF1(f_u, tau_n)
%EF1
% 1 - amostragem ideal
% 2 - Psi obtido pela fase do sinal analitico
% 3 - f_u[n] obtida pelo gradient de Psi
% 4 - considera f_u = f_i
% 5 
%   - retira 5% das amostras proximas das bordas
%   - considera f_r = media de f_u, aproximada pela mediana
%   - estima f_1 = mediana de f_u até tau_n
%   - estima f_2 = mediana de f_u a partir de tau_n
% 6 - estima r_i[n] pelo gradiente de f_i[n]

    NSamples = length(f_u);
    br = 0.05; %80/480; % fraction of samples to be ignored 
    brn = floor(br*NSamples); % number of samples to be ignored
    brmask1 = [(brn+1:tau_n)];
    brmask2 = [(tau_n+1:NSamples-brn)];
    brmask = [brmask1 brmask2];
    f_u_trunc = f_u(brmask);
    
    ri = zeros(NSamples,1);
    ri(brmask) = gradient(f_u_trunc);
   
    f1 = median(f_u(brmask1));
   f2 = median(f_u(brmask2));

   tau = tau_n/NSamples;
   f_r = tau*f1 + (1-tau)*f2;
    


