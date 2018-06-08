% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

D = open('Steps_30_05.mat')

dt = 2.0e-4;
SampleRate = 1/dt; %Hz
p = 1;
Fps = 10;  %frames per second
window = floor(SampleRate/Fps);
q = window;

p = 1 + 4*window + floor(window/2);
q = p + window - 1;
WholeSignal = [
D.Dadost1.MS_p_0';
D.Dadost2.MS_p_0';
D.Dadost3.MS_p_0';
D.Dadost4.MS_p_0';
D.Dadost5.MS_p_120';   %%%%%
D.Dadost6.MS_p_0';
D.Dadost7.MS_p_0';
D.Dadost8.MS_p_0';
D.Dadost9.PS_p_120';
];

%nominal values
F0 = 60; %Hz nominal
F1 = 60.0003; %Hz fundamental
AnalysisCycles = 6;
Ps = 120; %phase in degrees
Ph = Ps*pi/180;% Phase in radians
KaS = 0;
KxS = 0.1;
Vm = 1;

%Ph_corr_sec = 3.500; %us
%Ph_corr_deg = Ph_corr_sec*1e-6*F1;
Ph_corr_deg = 0.0;
Ph_corr = Ph_corr_deg*pi/180;
Mag_corr = 0.90563;
Vm = Vm*Mag_corr;
Ph = Ph + Ph_corr;

for k=1:5   %size(WholeSignal,1)           

    Signal = WholeSignal(k,p:q);
    NSamples = length(Signal);
    T = NSamples*dt;
    tau_nom = k*0.1*T;  %nominal tau
    
    %%%% Estimation of tau
    %Hilbert filter
    br = floor(0.05*NSamples); % 5% of NSamples are taken off at the beggining and end
    z=hilbert(Signal');  % calculates the analytic signal associated with Signal
    f_hi=unwrap(angle(z(2:end,:).*conj(z(1:end-1,:))));  % Hilbert estimate of the instantaneous frequency of Signal
    f=f_hi-median(f_hi(br:end-br));
    [ifmax, imax] = max(abs(f(br:end-br)));
    %threshold
    % th = 0.05;
    % imax_th = find(abs(f(br:end-br))>th,1);
    % imax_th = imax_th + br - 1;
    imax_th = imax + br - 1;
    if imax_th>0
        tau_e = imax_th*dt;
    else
        tau_e = 0;
        tau_nom = 0;
        imax_th = size(t,2)+1;
    end
    tau_error(k) = (tau_e - tau_nom);

    %Time vector using nominal tau
    tau_pp = tau_nom/T;
    n = -NSamples/2:(NSamples/2-1); %discrete time vector
    n = n - (tau_pp - 0.5)*NSamples;
    t = n*dt; %time vector
    u = zeros(1,length(t));
    tau = tau_e - tau_pp*T;
    u(1,t >= tau) = u(1,t >= tau) + 1;
    
    if KaS ~= 0  %phase step
        %f = @(x) x(1)*cos(x(2)*t + x(3) + x(4)*(pi/180)*u);
        f = @(x) x(1)*cos(x(2)*t + x(3) + x(4)*(pi/180)*u);
        xr = [Vm 2*pi*F1 Ph KaS];
    else    %mag step
        %f = @(x) x(1)*(1+x(2)*u).*cos(x(3)*t + x(4));
        f = @(x) x(1)*(1+x(2)*u).*cos(2*pi*F1*t + x(4));
        xr = [Vm KxS 2*pi*F1 Ph];
    end
    Ynom = f(xr);
    err = @(x) (Signal - f(x));
    
        % Levenberg-marquardt from optimization toolbox
    tol = 1e-7;
    OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','OptimalityTolerance',tol);
    OPTIONS.StepTolerance = 1e-12;
    [X,RESNORM,RESIDUAL,exitflag,output] = lsqnonlin(err,xr,[],[],OPTIONS);
    Y = f(X);

    plot(t,Signal,'b.',t,Y,'r',t,Ynom,'g'); legend('Signal','Y_{LM}','Y_{nom}')
    
    %Fasor medio
    T = NSamples*dt;
    tau_est = tau_pp + tau_error(k)/T;
    if KaS ~= 0
        %phase step
        Xe_r(k) = xr(1);
        Xe(k) = X(1);
        Phe_r(k) = xr(3)*tau_pp + (xr(3) + xr(4)*(pi/180))*(1 - tau_pp);  % [rad]
        Phe_r_deg(k) = Phe_r(k)*180/pi; % [deg]
        Phe(k) = X(3)*tau_est + (X(3)+X(4)*(pi/180))*(1 - tau_est);       % [rad]
        Freq(k) = X(2)/(2*pi);
    else
        %mag step
        Xe_r(k) = xr(1)*tau_pp + xr(1)*(1+xr(2))*(1 - tau_pp);
        Xe(k) = X(1)*tau_est + X(1)*(1+X(2))*(1-tau_est);
        Phe_r(k) = xr(4);   %  [rad]
        Phe_r_deg(k) = xr(4)*180/pi; % [deg]
        Phe(k) = X(4);   %  [rad]
        Freq(k) = X(3)/(2*pi);
    end

    if Ps>350
        Phe_deg = (Phe - 2*pi)*180/pi;
        Phe_r_deg(k) = Phe_r_deg(k) - 360;
    else
        Phe_deg = Phe*180/pi;
    end


    Phasor_mag_error(k,:) = (Xe(k) - Xe_r(k))/Xe_r(k);
    Phasor_ph_error(k,:) = (Phe(k) - Phe_r(k))/Phe_r(k);
    errors_perc(k,:) = 100*(xr - X)./xr;

    
end


