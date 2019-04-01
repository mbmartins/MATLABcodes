function [tau_error,extremos,det_mag,det_fase,det_nan] = MC_tau_error_hibrid_detector_v4(SNR,Pin,lim_mag,lim_fase,tau_pp_set,Nruns)
%Detection of tau using hilbert transform
display('v4 running')

UF1 = 0.05; %uncertainty of frequency in [%]
UVm = 1; %uncertainty of magnitude in [%]
UKa = 1; %uncertainty of Ka in [%]
UKx = 1; %uncertainty of Kx in [%]
UPs = 1; %uncertainty of phase in degrees
%Pin = 90; %initial phase in degrees  %worst case is 90 degrees

%SNR = 50;  %signal noise ratio in dB
rng('shuffle');

extremos = 0;
det_mag = 0;
det_fase = 0;
det_nan = 0;
%lim_mag = 0.0;
%lim_fase = 8;

%Monte carlo for tau_error
for k = 1:Nruns
    %signal generation
    F0 = 60; F1 = 60; %[Hz]
    F1 = F1 + F1*2*UF1/100*(rand-0.5);
    SampleRate = 5000; %Hz
    dt = 1/SampleRate;
    AnalysisCycles = 6;
    NSamples = floor(AnalysisCycles*SampleRate/F0);
    n = -NSamples/2:(NSamples/2-1); %discrete time vector
    %tau_pp_set = 0.1:0.1:0.9; % relative time of step in percent of total time 
    tau_pp = tau_pp_set(randi(size(tau_pp_set,2)));
    tau_0 = (tau_pp - 0.5)*NSamples; %discrete time displacement
    n = n - tau_0;
    t = n*dt; %time vector
    Vm = 0.9; 
    Vm = Vm + Vm*2*UVm/100*(rand-0.5);
    Ps = Pin + 2*UPs/100*(rand-0.5);%phase in degrees
    Ps_deg(k) = Ps;

    % Phase in radians
    Ph = Ps*pi/180;

    KaS = 0;   % IEEE Std phase (angle) step index: 10 degrees
    KaS = KaS + KaS*2*UKa/100*(rand-0.5);
    KxS = 0.1;   % magnitude step index: 0.1 
    KxS = KxS + KxS*2*UKx/100*(rand-0.5);
    Wf = 2*pi*F1;  % fundamental frequency

    Xm = Vm; %for now, single phase; TODO: 6-channels
    Ain = zeros(length(Xm),length(t));
    % Amplitude Step: applied after time passes 0
    i = 1;
    Ain(i,:) = Xm(i);
    Ain(i,t >= 0) = Ain(i,t >= 0) * (1 + KxS(i));
    %Phase step
    Theta(i,:) = (Wf(i)*t) ...                         % Fundamental
                     + Ph(i);               % phase shift
    Theta(i,t >= 0) = Theta(i,t >= 0) + (KaS(i) * pi/180);
    cSignal = (Ain.*exp(-1i.*Theta));
     %dB SNR = 20 log_10 Asinal/Aruido => Aruido = Asinal/10^(SNR/20)
    %Anoise = Vm/(sqrt(2))/10^(SNR/20);
    rSignal = real(cSignal);
    var_noise = ((std(rSignal))/(10^(SNR/20)))^2;
    std_noise = (std(rSignal))/(10^(SNR/20));
    noise = std_noise*randn(1,length(rSignal));  
    %noise_un = std_noise*rand(1,length(rSignal));

    Signal =  rSignal + noise;
    %SNR_hist(k) = snr(Signal,noise);
    SNR_est(k)=10*log10(var(rSignal)/var(noise));
    %plot(Signal)
    %%%% Estimation of tau
    %br = 0.05*NSamples; % 5% of NSamples are taken off at the beggining and end
    br = 0.02*NSamples; % 2% of NSamples are taken off at the beggining and end
    z=hilbert(Signal');  % calculates the analytic signal associated with Signal

    df=gradient(unwrap(angle(z)));% Hilbert estimate of the instantaneous frequency of z
    df=abs(df-median(df(br:end-br))); %v3mod
    
    Ain = (Ain - mean(Ain))./abs(Ain);
    nn = 1:(NSamples);
   
    [ifmax_fase, imax_fase] = max(abs(df(br:end-br))); %max indices
    limiar_fase=lim_fase*median(abs(df(br:end-br)));
    
    gmi=gradient(abs(z)); % estimate of the instantaneous magnitude of z
    gmi=abs(gmi-median(gmi(br:end-br))); %v3mod
    [ifmax_mag, imax_mag] = max(abs(gmi(br:end-br))); %max indices
    limiar_mag=lim_mag*median(abs(gmi(br:end-br)));
    
    % Threshold detection - v4mod
    crit1 = ifmax_fase(1)/limiar_fase; %fase
    crit2 = ifmax_mag(1)/limiar_mag; %mag
    if (crit1<1 && crit2<1) %not valid detections
        tau_e=NaN;
        det_nan = det_nan+1;
    elseif crit1>crit2  %fase detec
        tau_e=(br + imax_fase(1)-1)*dt;
        det_fase = det_fase +1;
    else %mag detec
        tau_e=(br + imax_mag(1)-1)*dt;
        det_mag = det_mag + 1;
    end

    tau = (tau_0 + NSamples/2)*dt;
    tau_error(k) = tau_e - tau;
    
    if (abs(tau_error(k)/dt) > 2.1 || isnan(tau_error(k)))
        extremos = extremos + 1;
%         figure
%         subplot(2,1,1)
%         plot(nn,abs(df));
%         xlabel('Samples', 'Fontsize',17);ylabel('$\mid{d[n]}\mid$','Interpreter','latex','Fontsize',17)
%         subplot(2,1,2)
%         plot(nn,Signal); xlabel('Samples', 'Fontsize',17); ylabel('y(t) [V]', 'Fontsize',17)
      
    end 
end


% 
% max_dt = max(tau_error)/dt
% min_dt = min(tau_error)/dt

%figure
%hist((tau_error/dt),200)
%hist((tau_error/dt),-250:250)
% hist((tau_error/dt),-10:10)
% xlabel('\tau error [\Deltat]')
% ylabel('Occurrences')


% figure
% subplot(2,1,1)
% plot(nn,abs(f)); title('Instantaneous frequency (abs)')
% subplot(2,1,2)
% plot(t,Signal); title('Signal with Magnitude Step')