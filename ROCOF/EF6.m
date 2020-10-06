function [f1,f2,F,f_i,ri] = EF6(f_i, az, tau_n)
% EF6
% 1 - amostragem ideal
% 2 - Psi obtido pela fase do sinal analitico
% 3 - f_i[n] obtida pelo gradient de Psi
% 4 - considera f_u = polinomio p obtido pela aplicacao de PATV (d=1,lambda?) a f_i 
% compensado pela amplitude az[n], excluindo amostras proximo as bordas e
% proximas a tau
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
    d = 1; lambda = 2.5; 
    f_icomp = f_i.*az./median(az);
    f_i2 = f_icomp(brn:end-brn-1);
    
    Nit = 20;
    [x, f_u, cost, u, v] = patv_MM(f_i2, d, lambda, Nit);
    f_i = f_u + x;  
    
    %exclui amostras proximas a tau
    % f_i ou f_u?
    f1 = median(f_i(1:(tau_n-brn)));
    f2 = median(f_i((tau_n-brn+1):end));
    ri = gradient(f_i); % ou usar diretamente u e v
    ru = gradient(f_u);
    %F = median([f1 f2]);
    F = median(f_u);
