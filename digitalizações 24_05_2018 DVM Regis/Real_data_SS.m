% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

D = open('Resultados_25_05.mat')

% p = 241 + 4*480 +1;
% q = 240 + 5*480 +1;

p = 1;
q = 480;

WholeSignal = [D.Dados.CM5';
];

% Nominal values
F0 = 60; %Hz nominal
F1 = 60; %Hz fundamental
Vm = 1; %70*sqrt(2) =~ 100;
Xm = Vm;
Ps = 360; %phase in degrees (use 360 instead of zero)
Ph = Ps*pi/180;% Phase in radians
Wf = 2*pi*F1;  % fundamental frequency
SampleRate = 4800; %Hz
dt = 1/SampleRate;

AnalysisCycles = 6;

for pin = 1:10
    p = (pin-1)*480+1
    q = p + 480 -1;
    Signal = WholeSignal(1,p:q);
    NSamples = length(Signal);
    T = NSamples*dt;
    tau_pp = 0; % relative time of step in percent of total time 
    tau_r = tau_pp*T;  % [s]
    tau_0 = (tau_pp - 0.5)*NSamples; %discrete time displacement
    n = -NSamples/2:(NSamples/2-1); %discrete time vector
    n = n - tau_0;
    t = n*dt; %time vector

    % Signal Model
    f = @(x) x(1)*cos(x(2)*t + x(3) + 2*pi);
    xnom = [Vm 2*pi*F1 Ph];

    err = @(x) (Signal - f(x));
    
    xr = xnom;
    % Levenberg-marquardt from optimization toolbox
        tol = 1e-7;
        OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','OptimalityTolerance',tol);
        OPTIONS.StepTolerance = 1e-12;
        [X,RESNORM,RESIDUAL,exitflag,output] = lsqnonlin(err,xr,[],[],OPTIONS);
        Y = f(X);
    
    %Steady State LS fit
    SignalParams(1) = F1;
    SignalParams(2) = 0;
    SignalParams(3) = 0;
    SignalParams(7) = 0;
    SignalParams(8) = 0;
    DelayCorr = 0;
    MagCorr = 1;
    [Synx(pin,:),Freq(pin,:),ROCOF(pin,:)] = SteadyStateFit ( ...
                        SignalParams, ...
                        DelayCorr, ...
                        MagCorr, ...
                        F0, ...
                        AnalysisCycles, ...
                        SampleRate, ...
                        Signal ...
                    )
   
    Freq_LM(pin) = X(2)/(2*pi);
    
    errors(pin,:) = (X - xr)*100./xr; % error in [%]
    
    Mag_SS(pin,:) = abs(Synx(pin))*sqrt(2);
    Mag_LM(pin,:) = X(1);
    
    Phase_SS(pin,:) = angle(Synx(pin))*180/pi;  %in degrees
    Phase_LM(pin,:) = (X(3)-2*pi)*180/pi; % in degrees
    
    Delay_SS(pin,:) = Phase_SS(pin)/(360*Freq(pin)) % in [s]
end
        



