% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
clear all; close all; clc;

D = open('Digitalizacoes_02_05_2018.mat')

p = 241 + 4*480;
q = 240 + 5*480;
WholeSignal = [D.Dig_t0.SS_0';
%                 D.Dig_t1.PS_n_0';
%                 D.Dig_t2.PS_n_0';
%                 D.Dig_t3.PS_n_0';
%                 D.Dig_t4.PS_n_0';
%                 D.Dig_t5.PS_n_0';
%                 D.Dig_t6.PS_n_0';
%                 D.Dig_t7.PS_n_0';
%                 D.Dig_t8.PS_n_0';
%                 D.Dig_t9.PS_n_0';
];

for i=1:size(WholeSignal,1)           

    Signal = WholeSignal(i,p:q);
    NSamples = length(Signal);
    n = 1:NSamples;
    SampleRate = 4800; %Hz
    dt = 1/SampleRate;
    t = n*dt;
    T = NSamples*dt;
    tau = i*0.1*T;
    F0 = 60; %Hz nominal
    F1 = 60; %Hz fundamental
    AnalysisCycles = 6;

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
        tau = 0;
        imax_th = size(t,2)+1;
    end

    tau_error(i) = (tau_e - tau);
    
    
    %Frequency estimation
    %splittig Psi in two
    Psi = unwrap(angle(z));
    Psi_1 = Psi(br:imax_th-1-br);
    t_1 = t(br:imax_th-1-br);
    Psi_2 = Psi(imax_th+br:end-br);
    t_2 = t(imax_th+br:end-br);

    %estimation of angular frequency (w = 2*pi*f) using hilbert
    % w is the slope of the linear curve

    %1st degree model matrix
    %calculating 2freqs
%     X1 = [t_1' ones(size(t_1))']; y = Psi_1;
%     beta = (X1'*X1)\X1'*y;
%     fest_1 = beta(1)/(2*pi)
%     ferr_1 = 100*(fest_1 - F1)/F1
%     res_1 = (X1*beta - y);
% 
%     X2 = [t_2' ones(size(t_2))']; y = Psi_2;
%     beta = (X2'*X2)\X2'*y;
%     fest_2 = beta(1)/(2*pi)
%     ferr_2 = 100*(fest_2 - F1)/F1
%     res_2 = (X2*beta - y);
% 
%     fest = mean([fest_1 fest_2])
%     ferr = 100*(fest - F1)/F1
% 
%     figure
%     subplot(3,1,1); plot(t_1,Psi_1,'.',t_2,Psi_2,'.'); title('\Psi')
%     subplot(3,1,2); plot(t_1,res_1,'.',t_2,res_2,'.'); title('residuals = (X*\beta - y)')
%     subplot(3,1,3); plot(t(2:end),abs(f_hi)); title('Instantaneous frequency')

    %Other approach
    %calculating 1freq
    Psi = [Psi_1; Psi_2]; y = Psi;
    t_m = [t_1 t_2]';
    X = [t_m'; ones(size(t_m))']';
    beta = (X'*X)\X'*y;
    fest_m = beta(1)/(2*pi);
    ferr_m(i,1) = 100*(fest_m - F1)/F1;
    res_m = (X*beta - y);
%    figure
%     subplot(3,1,1); plot(t_m,Psi); title('\Psi')
%     subplot(3,1,2); plot(t_m,res_m,'.'); title('residuals = (X*\beta - y)')
%     subplot(3,1,3); plot(t(2:end),abs(f_hi)); title('Instantaneous frequency')
    
end
% subplot(2,1,1)
% plot(t(2:end),abs(f)); title('Instantaneous frequency')
% subplot(2,1,2)
% plot(t,Signal); title('Signal')


