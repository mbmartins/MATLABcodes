function [f1,f2,F,f_u,ri] = EF4(Psi_noise, az, Fs, tau_n,lambda)    
% EF4
% 1 - amostragem ideal
% 2 - Psi obtido pela fase do sinal analitico
% 3 - f_i[n] obtida pelo gradient de Psi filtrado por PATV (d=1, lambda
% para minimizar o FE de fr medio, Psi_i = p+x (filtro somente do ruido)
% 4 - considera f_u = f_i sem os saltos, retirados pela mediana em (5)
% 5 
%   - considera f_r = media de f_u, aproximada pela mediana
%   - estima f_1 = mediana de f_u até tau_n
%   - estima f_2 = mediana de f_u a partir de tau_n
% 6 - estima r_i[n] pelo gradiente de f_i[n]

    br = 0.05;%80/480; % fraction of samples to be ignored 
    brn = floor(br*NSamples)+1; % number of samples to be ignored
    brmask = [brn+1:NSamples-brn];
    
    Nit = 20; 
   
    [x, p, cost, u, v] = patv_MM(Psi_noise, d, lambda, Nit);
    %[x, p, cost, u, v] = patv_MM(Psi_i(brmask), d, lambda, Nit); 
    % OBS: nao funciona bem retirar as amostras inciais e finais antes do 
    % PATV para este estimador
    
    Psi_u = p;
    Psi_i = p + x;
    f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    f_u=gradient(Psi_u)*Fs/(2*pi);% Hilbert estimate of the underlying frequency of z   

    ri = zeros(1,NSamples);
    ri(brmask) = gradient(f_i(brmask));
    ru = gradient(f_u);
    
     f1 = median(f_i(1:tau_n));
     f2 = median(f_i((tau_n+1):end));
%     F = median(f_u); 
%f1 = mean(f_i(1:tau_n));
%f2 = mean(f_i((tau_n+1):end));
F = mean(f_i);

    phi_0_est = Psi_u(1)*180/pi;
