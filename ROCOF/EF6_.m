function [f1,f2,F,f_u,ri] = EF6(f_u, az, tau_n)
    NSamples = length(f_u);
    br = 80/480; % fraction of samples to be ignored 
    brn = floor(br*NSamples); % number of samples to be ignored
    brmask = [(brn+1:tau_n-brn+1) (tau_n+brn+2:NSamples-brn)];

%EF5    
    %retira as amostras proximo as bordas
    % considera a mediana de az sem as bordas
    az5 = az(brmask)/median(az(brmask));
    f_u = f_u(brmask).*az5;
    
    %retira as amostras proximo ao tau
    f1 = median(f_u(1:(tau_n-2*brn)));
    f2 = median(f_u((tau_n-2*brn+1):end));
    ri = gradient(f_u);  % OBS: a função gradiente 'suaviza' o degrau
    F = median([f1 f2]);