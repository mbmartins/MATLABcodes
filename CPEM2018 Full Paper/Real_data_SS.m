% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

D = open('Steps_complete.mat')

% p = 241 + 4*480 +1;
% q = 240 + 5*480 +1;

WholeSignal = [ D.Dadost1.PS_p_120';
                D.Dadost2.PS_p_120';
                D.Dadost3.PS_p_120';
                D.Dadost4.PS_p_120';
                D.Dadost5.PS_p_120';
                D.Dadost6.PS_p_120';
                D.Dadost7.PS_p_120';
                D.Dadost8.PS_p_120';
                D.Dadost9.PS_p_120';
];

% Nominal values
F0 = 60; %Hz nominal
F1 = F0; %Hz fundamental
Vm = 1; %70*sqrt(2) =~ 100;
Xm = Vm;
%Ps = 360; %phase in degrees (use 360 instead of zero)
%Ph = Ps*pi/180;% Phase in radians
Wf = 2*pi*F1;  % fundamental frequency

%dt = 2.0833e-4;  %-> observar que o DMM trunca neste valor!!!!
dt = 2.0e-4;
SampleRate = 1/dt; %Hz

p = 1;
Fps = 10;  %frames per second
window = floor(SampleRate/Fps);
q = window;

AnalysisCycles = 6;
for s = 1:9
for pin = 1:Fps-1
    p = (pin-1)*window+window/2+1;
    q = p + window -1;
    Signal = WholeSignal(s,p:q);
    NSamples = length(Signal);
%     T = NSamples*dt;
%     tau_pp = 0; % relative time of step in percent of total time 
%     tau_r = tau_pp*T;  % [s]
%     tau_0 = (tau_pp - 0.5)*NSamples; %discrete time displacement
%     n = -NSamples/2:(NSamples/2-1); %discrete time vector
%     n = n - tau_0;
%     t = n*dt; %time vector
% 
%     % Signal Model
%     f = @(x) x(1)*cos(x(2)*t + x(3) + 2*pi);
%     xnom = [Vm 2*pi*F1 Ph];
% 
%     err = @(x) (Signal - f(x));
%     
%     xr = xnom;
    % Levenberg-marquardt from optimization toolbox
%         tol = 1e-7;
%         OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','OptimalityTolerance',tol);
%         OPTIONS.StepTolerance = 1e-12;
%        [X,RESNORM,RESIDUAL,exitflag,output] = lsqnonlin(err,xr,[],[],OPTIONS);
%        Y = f(X);
    
    %Steady State LS fit
    SignalParams(1) = F1; %F1;
    SignalParams(2) = 0;
    SignalParams(3) = 0;
    SignalParams(7) = 0;
    SignalParams(8) = 0;
    DelayCorr = 0; %-3500.0; %   %[in nanosecond]
    MagCorr = 1;
    [Synx(pin,:),Freq(pin,s),ROCOF(pin,:)] = SteadyStateFit ( ...
                        SignalParams, ...
                        DelayCorr, ...
                        MagCorr, ...
                        F0, ...
                        AnalysisCycles, ...
                        SampleRate, ...
                        Signal ...
                    )
   
%    Freq_LM(pin) = X(2)/(2*pi);
    
%    errors(pin,:) = (X - xr)*100./xr; % error in [%]
    
    Mag_SS(pin,s) = abs(Synx(pin))*sqrt(2);
%    Mag_LM(pin,:) = X(1);
    
    Phase_SS(pin,s) = angle(Synx(pin))*180/pi;  %in degrees
%    Phase_LM(pin,:) = (X(3)-2*pi)*180/pi; % in degrees
    
    Delay_SS(pin,:) = Phase_SS(pin)/(360*Freq(pin)) % in [s]
    %Delay_LM(pin,:) = Phase_LM(pin)/(360*Freq_LM(pin)) % in [s]
end
end

mean_delay = mean(Delay_SS)

delay_corr = Delay_SS - mean(Delay_SS)

% tn = 0:dt:NSamples*dt-dt;
% Signal_est = Mag_SS(pin,:)*cos(2*pi*Freq(pin,:).*tn + Phase_SS(pin,:)*pi/180)
% Signal_nom = Vm*cos(2*pi*60.*tn + Ph)
% plot(tn,Signal_est-Signal)



