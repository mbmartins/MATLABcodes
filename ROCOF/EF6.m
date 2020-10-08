function [f1,f2,F,f_u,ri] = EF6(f_i, az, tau_n,lambda)
% EF6
% 1 - amostragem ideal
% 2 - Psi obtido pela fase do sinal analitico
% 3 - f_i[n] obtida pelo gradient de Psi
% 4 - considera f_u = polinomio p obtido pela aplicacao de PATV (d=1,lambda?) a f_i 
% compensado pela amplitude az[n],{EM REVISAO excluindo amostras proximo as bordas e
% proximas a tau EM REVISAO }
% 5 
%   - considera f_r = media de f_u, aproximada pela mediana
%   - estima f_1 = mediana de f_i até tau_n
%   - estima f_2 = mediana de f_i a partir de tau_n
% 6 - estima r_i[n] pelo gradiente de f_i[n]
%   - estima r_u[n] pelo gradiente de f_u[n]

NSamples = length(f_i);
    br = 0.05; %80/480; % fraction of samples to be ignored 
    brn = floor(br*NSamples); % number of samples to be ignored

    %Ef6 - PATV aplicado a fi
    d = 1; %lambda = 2.5; 
    f_icomp = f_i.*az./median(az);
    f_i2 = f_icomp(brn:end-brn-1);
    
    Nit = 20;
    [x, f_u, cost, u, v] = patv_MM(f_i2, d, lambda, Nit);
    f_i = f_u + x;  
    ri = gradient(f_i); % ou usar diretamente u e v
    ru = gradient(f_u);

    
    
%      f1 = median(f_i(1:(tau_n)));
%      f2 = median(f_i((tau_n+1):end));
  %F = median(f_u);

      f1 = mean(f_i(1:(tau_n)));
      f2 = mean(f_i((tau_n+1):end));
  
    F = mean(f_u);