function [freq] = MedFR(Signal,Fs)
% estimates freq [Hz] from the input Signal
% based on Hilberts analytic signal
NSamples = size(Signal,1);
    br = floor(0.05*NSamples); % 5% of NSamples are taken off at the beggining and end
    z=hilbert(Signal');  % calculates the analytic signal associated with Signal
    df=gradient(unwrap(angle(z)));% Hilbert estimate of the instantaneous frequency of z
    fest = median(df(br:end-br)); %frequency estimate
    freq = fest*Fs/(2*pi); % in Hz
   