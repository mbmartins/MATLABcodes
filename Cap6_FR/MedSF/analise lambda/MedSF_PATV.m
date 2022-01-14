function [f1,f2,f_r,f_u,ri,ru] = MedSF_PATV(f_i, az, tau_n, lambda)
% MedSF_PATV

    NSamples = length(f_i);
    br = 0.05;%80/480; % fraction of samples to be ignored 
    q = floor(br*NSamples); % number of samples to be ignored
    brmask = [q+1:NSamples-q];
    
    d = 0; 
    
    %exclui as amostras proximo as bordas
    f_i2 = f_i(brmask);
    
    Nit = 200;
    [x, f_u, cost, u, v] = patv_MM(f_i2, d, lambda, Nit);
    f_iest = f_u + x; 
    ri = zeros(NSamples,1);ru = zeros(NSamples,1);
    ri(brmask) = gradient(f_iest); % ou usar diretamente u
    
    % indices deslocados pela supressao das amostras nas bordas
    % tratamento de erro para o caso de tau ser muito proximo a borda
    n1 = (tau_n-q);
    if n1<1
        n1 = 1;
    end
    n2 = (tau_n-q+1);
    if n2<1
        n2 = 1;
    end;
    
   f1 = median(f_iest(1:n1));
   f2 = median(f_iest(n2:end));

      tau = tau_n/NSamples;
    f_r = tau*f1 + (1-tau)*f2;
