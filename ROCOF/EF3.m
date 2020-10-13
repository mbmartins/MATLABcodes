function [f1,f2,F,f_u,ri] = EF3(f_i, az, tau_n)
%EF3
% 1 - amostragem ideal
% 2 - Psi obtido pela fase do sinal analitico
% 3 - f_i[n] obtida pelo gradient de Psi
% 4 - considera f_u = f_i compensado pela amplitude az[n], ja com as
% amostras proximas das bordas retiradas
% 5 
%   - retira 5% das amostras proximas das bordas
%   - retira 5% das amostras proximas de tau
%   - considera f_r = media de f_u, aproximada pela mediana
%   - estima f_1 = mediana de f_u até tau_n
%   - estima f_2 = mediana de f_u a partir de tau_n
% 6 - estima r_i[n] pelo gradiente de f_i[n]

    NSamples = length(f_i);
    br = 0.05; %80/480; % fraction of samples to be ignored 
    brn = floor(br*NSamples); % number of samples to be ignored
    brmask1 = [(brn+1:tau_n-brn)];
    brmask2 = [(tau_n+brn+1:NSamples-brn)];
    brmask = [brmask1 brmask2];

    %retira amostras próximas a borda
    aztrunc = az(brmask);
    f_utrunc = f_i(brmask);
    
    f_u = f_utrunc.*aztrunc./median(aztrunc);
    ri = gradient(f_i);

    %retirando amostras próximas a tau
    % 2*brn pois no f_u2 já foi retirado brn
%      f1 = median(f_u(1:(tau_n-2*brn)));  
%      f2 = median(f_u((tau_n-2*brn+1):end));
%     F = median(f_u);
  f1 = mean(f_u(1:(tau_n-2*brn)));
  if (tau_n-2*brn+1)<length(f_u)
    f2 = mean(f_u((tau_n-2*brn+1):end));
  else
    f2 = f_u(end);
  end
  F = mean(f_u);