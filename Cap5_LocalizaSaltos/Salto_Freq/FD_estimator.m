function [tau_e,dmax,limiar] = FD_estimator(Signal,kf)
% estima a localiza��o do salto de frequ�ncia do sinal

N = length(Signal);
br = floor(0.1*N); % 5% of NSamples are taken off at the beggining and end
%excluindo amostras do inicio e do final
brmask = br:N-br-1;

z=hilbert(Signal');  % calculates the analytic signal associated with Signal
fi=gradient(unwrap(angle(z)));% Hilbert estimate of the instantaneous frequency of z
ai = abs(z);
ai = ai(brmask)/median(ai); %normaliza��o
fi = fi(brmask);
%compensa��o pela magnitude
%fic = fi.*ai;

ri = gradient(fi); %rocof de fi compensado

d_r = abs(ri);

[dmax, imax_freq] = max(d_r); %max indices
limiar=kf*median(d_r); %limiar de detec��o em fun��o da mediana do sinal de detec��o
    
    % Threshold detection 
    if dmax > limiar 
        tau_e=(br + imax_freq(1)-1); %em delta t
    else
        tau_e = NaN;
    end

