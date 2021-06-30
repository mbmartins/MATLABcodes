function [tau,freq] = HE_tau(km,kf,Signal,F0,Fs)
% estimates tau [s] and freq [Hz] from the Signal
% based on Hilberts analytic signal
NSamples = size(Signal,1);
dt = 1/Fs;
    %%%% Estimation of tau
    br = 0.05*NSamples; % 5% of NSamples are taken off at the beggining and end
    %br = 0.02*NSamples; % 2% of NSamples are taken off at the beggining and end
    z=hilbert(Signal');  % calculates the analytic signal associated with Signal

    df=gradient(unwrap(angle(z)));% Hilbert estimate of the instantaneous frequency of z
    fest = median(df(br:end-br)); %frequency estimate
    freq = fest*Fs/(2*pi);
    df=abs(df-fest); %v3mod
    
    %Ain = (Ain - mean(Ain))./abs(Ain);
    nn = 1:(NSamples);
   
    [ifmax_fase, imax_fase] = max(abs(df(br:end-br))); %max indices
    limiar_fase=kf*median(abs(df(br:end-br)));
    
    gmi=gradient(abs(z)); % estimate of the gradient of instantaneous magnitude of z
    gmi=abs(gmi-median(gmi(br:end-br))); %v3mod
    [ifmax_mag, imax_mag] = max(abs(gmi(br:end-br))); %max indices
    limiar_mag=km*median(abs(gmi(br:end-br)));
    
    % Threshold detection 
    crit1 = ifmax_fase(1)/limiar_fase; %fase
    crit2 = ifmax_mag(1)/limiar_mag; %mag
    if (crit1<1 && crit2<1) %not valid detections
        tau=NaN;
        %det_nan = det_nan+1;
    elseif crit1>crit2  %fase detec
        tau=(br + imax_fase(1)-1)*dt;
        %det_fase = det_fase +1;
    else %mag detec
        tau=(br + imax_mag(1)-1)*dt;
        %det_mag = det_mag + 1;
    end