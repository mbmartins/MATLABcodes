function [tau_e,dmax,limiar] = FD_PATV_estimator(Signal,kf,lambda)
% estima a localização do salto de frequência do sinal
%fazer a estimacao de fi igual ao antigo EF6

N = length(Signal);
br = floor(0.05*N); % 5% of NSamples are taken off at the beggining and end

z=hilbert(Signal');  % calculates the analytic signal associated with Signal

fi=gradient(unwrap(angle(z)));% Hilbert estimate of the instantaneous frequency of z
%excluindo amostras do inicio e do final
brmask = br:N-br-1;
fi = fi(brmask);

%aplicar PATV em fi
d = 0; Nit = 100;
[x, p, cost, u, v] = patv_MM(fi, d, lambda, Nit);

ri = gradient(x); %rocof de x (a componente TV de fi)

d_r = abs(ri);
%plot(d_r);

[dmax, imax_freq] = max(d_r); %max indices
limiar=kf; %limiar de detecção
    
    % Threshold detection
    if dmax(1) > limiar 
        tau_e=(br + imax_freq(1)-1); %em [delta t]
    else
        tau_e = NaN;
    end

