% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

D = open('Digitalizacoes_02_05_2018.mat')

p = 241 + 4*480 +1;
q = 240 + 5*480 +1;

WholeSignal = [%D.Dig_t0.SS_0';
                 D.Dig_t1.MS_p_120';
                 D.Dig_t2.MS_p_120';
                 D.Dig_t3.MS_p_120';
                 D.Dig_t4.MS_p_120';
                 D.Dig_t5.MS_p_120';
                 D.Dig_t6.MS_p_120';
                 D.Dig_t7.MS_p_120';
                 D.Dig_t8.MS_p_120';
                 D.Dig_t9.MS_p_120';
];

% Nominal values
F0 = 60; %Hz nominal
F1 = 60; %Hz fundamental
Vm = 1; %70*sqrt(2) =~ 100;
Xm = Vm;
Ps = 120; %phase in degrees (use 360 instead of zero)
Ph = Ps*pi/180;% Phase in radians
KxS = 0.1;   % magnitude step index: 0.1 
Vm = Vm/(1 + KxS);
KaS = 0;   % IEEE Std phase (angle) step index: 10 degrees
Wf = 2*pi*F1;  % fundamental frequency
SampleRate = 4800; %Hz
dt = 1/SampleRate;

AnalysisCycles = 6;
Signal = WholeSignal(1,p:q);
NSamples = length(Signal);
T = NSamples*dt;


for ti=1:size(WholeSignal,1)           

    Signal = WholeSignal(ti,p:q);
    tau_pp = 0.1*ti; % relative time of step in percent of total time 
    tau_r = tau_pp*T;  % [s]
    tau_0 = (tau_pp - 0.5)*NSamples; %discrete time displacement
    n = -NSamples/2:(NSamples/2-1); %discrete time vector
    n = n - tau_0;
    t = n*dt; %time vector
    

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
            tau_e = imax_th*dt;   %[s]
        else
            tau_e = 0;
            tau_r = 0;
            imax_th = size(t,2)+1;
        end
        tau_error(ti) = (tau_e - tau_r);
    
    %%%%% Estimation of parameters
    
    u = zeros(length(Xm),length(t));
    u(length(Xm),t >= 0) = u(length(Xm),t >= 0) + 1;
    
    % Signal Model
    if KaS ~= 0  % Phase step % Model: f(x) = x1*cos(w*t + phi + x2*u(t - tau)) 
        f = @(x) x(1)*cos(x(2)*t + x(3) + 2*pi + x(4)*(pi/180)*u);
        xnom = [Vm 2*pi*F1 Ph KaS];
    else         % Mag Step    %f = @(x) x(1)*(1+x(4)*u).*cos(x(2)*t + x(3));
        f = @(x) x(1)*(1+x(2)*u).*cos(x(3)*t + x(4));
        xnom = [Vm KxS 2*pi*F1 Ph];
    end
    err = @(x) (Signal - f(x));
    
    xr = xnom;
    % Levenberg-marquardt from optimization toolbox
        tol = 1e-7;
        OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','OptimalityTolerance',tol);
        OPTIONS.StepTolerance = 1e-12;
        [X,RESNORM,RESIDUAL,exitflag,output] = lsqnonlin(err,xr,[],[],OPTIONS);
        Y = f(X);
     
    % Weighted mean reference Phasor
    tau_est = tau_e/T;
        if KaS ~= 0
            %phase step
            Xe_r = xr(1);
            Xe = X(1);
            Phe_r = xr(3)*tau_pp + (xr(3) + xr(4)*(pi/180))*(1 - tau_pp);
            Phe = X(3)*tau_est + (X(3)+X(4)*(pi/180))*(1 - tau_est);
        else
            %mag step
            Xe_r = xr(1)*tau_pp + xr(1)*(1+xr(2))*(1 - tau_pp);
            Xe = X(1)*tau_est + X(1)*(1+X(2))*(1-tau_est);
            Phe_r = xr(4);
            Phe = X(4);
        end
        
        Phasor_mag_error(ti,:) = (Xe - Xe_r)*100/Xe_r;  % error in [%]
        Phasor_ph_error(ti,:) = (Phe - Phe_r)*100/Phe_r; % error in [%]
        errors(ti,:) = (X - xr)*100./xr; % error in [%]
        
    
end



