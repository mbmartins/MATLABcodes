function [tau_e,dmax,limiar] = HF_estimator(Signal,kf)
% estima a localiza��o do salto de frequ�ncia do sinal

N = length(Signal);
br = floor(0.05*N); % 5% of NSamples are taken off at the beggining and end

z=hilbert(Signal');  % calculates the analytic signal associated with Signal

df=gradient(unwrap(angle(z)));% Hilbert estimate of the instantaneous frequency of z
ri = gradient(df); %rocof de fi

%excluindo amostras do inicio e do final
brmask = br:N-br-1;
df = df(brmask);
d_r = abs(ri(brmask));

[dmax, imax_freq] = max(d_r); %max indices
limiar=kf*median(d_r); %limiar de detec��o em fun��o da mediana
    
    % Threshold detection 
    if dmax > limiar 
        tau_e=(br + imax_freq(1)-1); %em delta t
    else
        tau_e = NaN;
    end

