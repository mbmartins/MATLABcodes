function [f1,f2,F,f_u,ri] = EF5(f_i, az, tau_n, lambda)
% EF5
% 1 - amostragem ideal
% 2 - Psi obtido pela fase do sinal analitico
% 3 - f_i[n] obtida pelo gradient de Psi
% 4 - considera f_u = polinomio p obtido pela aplicacao de PATV (d=1,lambda?) a f_i 
% compensado pela amplitude az[n], excluindo amostras proximo as bordas
% 5 
%   - considera f_r = media de f_u, aproximada pela mediana
%   - estima f_1 = mediana de f_i até tau_n
%   - estima f_2 = mediana de f_i a partir de tau_n
% 6 - estima r_i[n] pelo gradiente de f_i[n]
%   - estima r_u[n] pelo gradiente de f_u[n]

    NSamples = length(f_i);
    br = 0.05;%80/480; % fraction of samples to be ignored 
    brn = floor(br*NSamples)+1; % number of samples to be ignored

    d = 0; %lambda = 2.5; 
    % com d = 1, fu pode assumir variacao lenta
    f_icomp = f_i.*az./median(az);
    
    %exclui as amostras proximo as bordas
    f_i2 = f_icomp(brn:end-brn-1);
    %f_u2 = f_ucomp;
    
    Nit = 20;
    [x, f_u, cost, u, v] = patv_MM(f_i2, d, lambda, Nit);
    f_i = f_u + x; 
    ri = gradient(f_i); % ou usar diretamente u
    ru = gradient(f_u);
    
    % indices deslocados pela supressao das amostras nas bordas
    % tratamento de erro para o caso de tau ser muito proximo a borda
    n1 = (tau_n-brn);
    if n1<1
        n1 = 1;
    end
    n2 = (tau_n-brn+1);
    if n2<1
        n2 = 1;
    end;
    
     f1 = median(f_i(1:n1));
     f2 = median(f_i(n2:end));
%     F = median(f_i);
%     f1 = mean(f_i(1:n1));
%     f2 = mean(f_i(n2:end));

%    F = mean(f_i);
      tau = tau_n/NSamples;
    F = tau*f1 + (1-tau)*f2;