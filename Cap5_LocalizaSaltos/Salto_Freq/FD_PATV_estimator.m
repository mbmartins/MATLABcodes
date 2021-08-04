function [tau_e,dmax,Lr] = FD_PATV_estimator(Signal,Lr,lambda)
% estima a localização do salto de frequência do sinal

N = length(Signal);
br = floor(0.05*N); % 5% of NSamples are taken off at the beggining and end

z=hilbert(Signal');  % calculates the analytic signal associated with Signal

fi=gradient(unwrap(angle(z)));% Hilbert estimate of the instantaneous frequency of z
ai = abs(z);

%excluindo amostras do inicio e do final
brmask = br:N-br-1;
ai = ai(brmask)/median(ai); %normalização
fi = fi(brmask);
%compensação pela magnitude
fic = fi.*ai;

%aplicar PATV em fi
d = 0; Nit = 100;
[x, p, cost, u, v] = patv_MM(fi, d, lambda, Nit);

ri = diff(x); %rocof de x (a componente TV de fi)

% plot(fic/max(abs(fic))); 
% hold on; plot(ri/max(abs(ri))); legend('fic','ri')

d_r = abs(ri);
%plot(d_r);

[dmax, imax_freq] = max(d_r); %max indices
    
    % Threshold detection
    if dmax(1) > Lr
        tau_e=(br + imax_freq(1)-1); %em [delta t]
    else
        tau_e = NaN;
    end
    
%     %DEBUG
%     if tau_e - 240> 3
%         message = "erro = "+tau_e
%     end
    

