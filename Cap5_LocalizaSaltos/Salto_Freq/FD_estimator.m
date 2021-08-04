function [tau_e,dmax,limiar] = FD_estimator(Signal,kf)
% estima a localização do salto de frequência do sinal

N = length(Signal);
br = floor(0.1*N); % 5% of NSamples are taken off at the beggining and end
%excluindo amostras do inicio e do final
brmask = br:N-br-1;

z=hilbert(Signal');  % calculates the analytic signal associated with Signal
fi=gradient(unwrap(angle(z)));% Hilbert estimate of the instantaneous frequency of z
ai = abs(z);
ai = ai(brmask)/median(ai); %normalização
fi = fi(brmask);
%compensação pela magnitude
%fic = fi.*ai;

ri = gradient(fi); %rocof de fi compensado

d_r = abs(ri);

[dmax, imax_freq] = max(d_r); %max indices
limiar=kf*median(d_r); %limiar de detecção em função da mediana do sinal de detecção
    
    % Threshold detection 
    if dmax > limiar 
        tau_e=(br + imax_freq(1)-1); %em delta t
    else
        tau_e = NaN;
    end

