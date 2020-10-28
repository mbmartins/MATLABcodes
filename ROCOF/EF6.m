function [f1,f2,f_r,f_u,ri,ru] = EF6(f_i, az, tau_n,lambda)
% EF6
% 1 - amostragem ideal
% 2 - Psi obtido pela fase do sinal analitico
% 3 - f_i[n] obtida pelo gradient de Psi
% 4 - considera f_u = polinomio p obtido pela aplicacao de PATV (d=1,lambda?) a f_i 
% compensado pela amplitude az[n],{EM REVISAO excluindo amostras proximo as bordas e
% proximas a tau EM REVISAO - a exclusao de amostras quando d=1 nao é uma boa}
% 5 
%   - considera f_r = media de f_u, aproximada pela mediana
%   - estima f_1 = mediana de f_i até tau_n
%   - estima f_2 = mediana de f_i a partir de tau_n
% 6 - estima r_i[n] pelo gradiente de f_i[n]
%   - estima r_u[n] pelo gradiente de f_u[n]

NSamples = length(f_i);
    br = 0.05; %80/480; % fraction of samples to be ignored 
    brn = floor(br*NSamples); % number of samples to be ignored
    brmask = [brn+1:NSamples-brn];
    
    %Ef6 - PATV aplicado a fi
    d = 1; %lambda = 2.5; 
    %f_icomp = f_i.*az./median(az);
    % TESTE - DEBUG
    f_icomp = f_i;
    
    %f_i2 = f_icomp(brn:end-brn-1);
    f_i2 = f_icomp;
    
    Nit = 20;
    [x, f_u, cost, u, v] = patv_MM(f_i2, d, lambda, Nit);
    f_i = f_u + x;
    ri = zeros(NSamples,1);ru = zeros(NSamples,1);
    ri(brmask) = gradient(f_i(brmask)); % ou usar diretamente u e v
    ru(brmask) = gradient(f_u(brmask));

    % ---- Deteccao de tau
    d = abs(ri);
    th = 0.1;
    [dmax,tau_n_est] = max(d);

    %----- Debug -------
%     figure;
%     nn = brn:NSamples-brn-1;
%     plot(nn,ri);
%     xlabel('x','Interpreter','latex');
%     ylabel('y','Interpreter','latex');
%     xlabel('Samples')
%     ylabel('$\hat{r}_i[n]$ [Hz/s]')
    

      f1 = median(f_i(1:tau_n));
      %f1 = mean(f_i(1:(tau_n)));
      
      if tau_n-brn<length(f_i)
        f2 = median(f_i((tau_n+1):end));
        %f2 = mean(f_i((tau_n+1):end));
      else
        f2 = f_i(end);
      end
       %F = median(f_i);
       
      tau = tau_n/NSamples;
    f_r = tau*f1 + (1-tau)*f2;
