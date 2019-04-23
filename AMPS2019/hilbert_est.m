function [fi, gmi] = hilbert_est(Signal)

%    br = 0.02*NSamples; % 2% of NSamples are taken off at the beggining and end
    z=hilbert(Signal');  % calculates the analytic signal associated with Signal
    fi = unwrap(angle(z));
    gmi = abs(z);
    
    
    %df=gradient(fi);% Hilbert estimate of the instantaneous frequency of z
    %df=abs(df-median(df(br:end-br))); %v3mod
    
    %Ain = (Ain - mean(Ain))./abs(Ain);
    %nn = 1:(NSamples);
   
    %[ifmax_fase, imax_fase] = max(abs(df(br:end-br))); %max indices
    %limiar_fase=lim_fase*median(abs(df(br:end-br)));
    
    %gmi=gradient(abs(z)); % estimate of the instantaneous magnitude of z
    %gmi=abs(gmi-median(gmi(br:end-br))); %v3mod
    %[ifmax_mag, imax_mag] = max(abs(gmi(br:end-br))); %max indices
    %limiar_mag=lim_mag*median(abs(gmi(br:end-br)));