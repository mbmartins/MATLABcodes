function [FE_LM, FE_SS, TVE_LM, TVE_SS] = LM_comparison_FE_TVE(SNR,F1,Ps, hm, ha)

% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 
%clear all; close all; clc;

%Parameters for signal generation
F0 = 60; %Hz
%F1 = 60.1; %Hz
SampleRate = 5000; %Hz
dt = 1/SampleRate;
AnalysisCycles = 48;
NSamples = floor(AnalysisCycles*SampleRate/F0);
n = -NSamples/2:(NSamples/2-1); %discrete time vector

%Nominal parameters, change if you want
Vm = 1; %70*sqrt(2) =~ 100;
Xm = Vm;
%Ps = 120; %phase in degrees
Ph = Ps*pi/180;% Phase in radians
%KaS = 0;   % IEEE Std phase (angle) step index: 10 degree
% OBS: KaS is used to track whether the test is phase step (KaS ~= 0) or not (KaS == 0) 
%KxS = -0.1;   % magnitude step index: 0.1
Wf = 2*pi*F1;  % fundamental frequency
%SNR = 80; %dB SNR = 20 log_10 Asinal/Aruido => Aruido = Asinal/10^(SNR/20)

%uncertainties of parameters in signal generation
if ha ~= 0
% phase         X1  w    ph   x3 (KaS)
    par_var = [1   0.05 1 1]; % parameter variation in percent related to nominal
else
% mag          x1  x2(KxS)  wf    ph  
    par_var = [1   1     0.05  1]; % parameter variation in percent related to nominal            
end





%%%% parameters for the LM window
        LM_cycles = 12;
        LMi = (NSamples - LM_cycles*SampleRate/F0)/2 + 1;
        LMf = LMi + LM_cycles*SampleRate/F0;

ti = 5;     % change ti to control the tau parameter: 
            % ex: ti=1 -> tau = 10% of T
            % ex: ti=5 -> tau = 50% of T
            % OBS: the tau parameter is not estimated in this code, you
            % need to give a value for it

    tau_pp = 0.1*ti; % relative time of step in percent of total time 
    tau_0 = (tau_pp - 0.5)*NSamples; %discrete time displacement
    n = -NSamples/2:(NSamples/2-1); %discrete time vector
    n = n - tau_0;
    t = n*dt; %relative time vector
    % note that the relative time vector is displaced by tau

    %%%% Signal parameters estimation
    % Here, u represents a step function
    u = zeros(length(Xm),length(t));
    tau = 0;  % the step happens always when the relative time t=0
    u(length(Xm),t >= tau) = u(length(Xm),t >= tau) + 1; 
    
 %%%%%  Defining the Mathematical model functions: f(x)
    % f is a function that represents the mathematical model
    % we infer if the model is phase step from KaS parameter
    % x is a vector that contains the parameters, depending on the case
    % x = [x1,w,phi,x3] for phase step
    % x = [x1,x2,w,phi] for mag step
    % xnom contains the nominal values set in the beggining
    if ha ~= 0
        % Phase step
        % Model: f(x) = x1*cos(w*t + phi + x3*u(t - tau)) 
        f = @(x,lmi,lmf) x(1)*cos(x(2)*t(lmi:lmf) + x(3) + 2*pi + x(4)*(pi/180)*u(lmi:lmf));
        xnom = [Vm 2*pi*F1 Ph ha];
    else %if not, it is a magnitude step
        % Mag step
        % model: f(x) = x1(1+x2*u(t-tau))*cos(w*t + phi)
        f = @(x,lmi,lmf) x(1)*(1+x(2)*u(lmi:lmf)).*cos(x(3)*t(lmi:lmf) + x(4));
        xnom = [Vm hm 2*pi*F1 Ph];
    end
    
%%%%%% Generating the Signal %%%%%%%
    % xr is the vector of parameters used to generate the signal, the
    % 'actual' values that the fitting should reach
    % xr can be the nominal values or something slightly different, to
    % simulate the uncertainty in the generation
    xr = xnom;  %xr is fixed - first guess are the nominal values: change if you want
    rng('shuffle');
    rn = (rand(1,length(xnom))-0.5);
    xr = xnom.*(1+2*(par_var/100).*rn);
    k = 1;  %it would be the kth Monte Carlo run
 
    %%%%%% Signal parameters %%%%%
            p = 1; Xm = xr(1); 
            if ha ~= 0
                Wf(p) = xr(2); 
                Ph(p) = xr(3);
                KxS_(p) = 0;
                KaS_(p) = xr(4);
            else
                Ph(p) = xr(4); 
                Wf(p) = xr(3);
                KxS_(p) = xr(2);
                KaS_(p) = 0;
            end
            
       %%%%% Signal generation
            Ain = zeros(length(Xm),length(t));
            % Amplitude Step: applied after time passes 0
            Ain(p,:) = Xm(p);
            Ain(p,t >= 0) = Ain(p,t >= 0) * (1 + KxS_(p));
            Wfv(1:length(t)) = Wf;
            %Wfv(1251:1500) = Wf + 0.008*2*pi;
            %Wfv(2251:2501) = Wf + 0.004*2*pi;
            
            %Phase step
            Theta(p,:) = (Wfv.*t) ...                         % Fundamental
                             + Ph(p);               % phase shift
            Theta(p,t >= 0) = Theta(p,t >= 0) + (KaS_(p) * pi/180);
            cSignal = (Ain.*exp(-1i.*Theta));
            rSignal = real(cSignal);
            
            % noise generation, given a SNR
            var_noise = ((std(rSignal))/(10^(SNR/20)))^2;
            std_noise = (std(rSignal))/(10^(SNR/20));
            noise = std_noise*randn(1,length(rSignal));
            
            Signal = rSignal + noise;
%            SNR_hist = snr(Signal,noise);  %just to check the SNR, if you want
            
         %y0: signal generated with the nominal parameters, without noise
         y0 = f(xnom,LMi,LMf);
%          plot(n(LMi:LMf),Signal(LMi:LMf),'.b',n(LMi:LMf),y0,'r')
%          title('Signal for the LM estimation window')
%          legend('Sampled Signal', 'Signal Model')
         
         % in this plot, if we used the nominal parameters to generate Signal, what we see is just the noise
%          figure
%          plot(n,Signal - y0,'r')
%          title('Point by point estimation error: Signal - y')

  %%%%%% Non linear fitting: now we finally get the estimated parameters
        % Levenberg-marquardt from optimization toolbox
        % Define the error function will to be used in the parameter estimation
        % OBS: you don't need to make it err^2 here

        % this applies only for the centered window
%        err2 = @(x) (Signal(LMi-250:LMf-251) - f(x,LMi-250,LMf-251));
        err5 = @(x) (Signal(LMi:LMf) - f(x,LMi,LMf));
%        err7 = @(x) (Signal(LMf-250:LMf+249) - f(x,LMf-250,LMf+249));
        x0 = xnom; %x0 is the first guess for the NL fitting
        tol = 1e-7;
        OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','OptimalityTolerance',tol);
        OPTIONS.StepTolerance = 1e-12;
                    OPTIONS.Display = 'off';
%        [X2,RESNORM,RESIDUAL,exitflag,output] = lsqnonlin(err2,x0,[],[],OPTIONS);
        [X5,RESNORM,RESIDUAL,exitflag,output] = lsqnonlin(err5,x0,[],[],OPTIONS);
%        [X7,RESNORM,RESIDUAL,exitflag,output] = lsqnonlin(err7,x0,[],[],OPTIONS);
        % X receives the estimated parameters
        % Y is the signal generated with the estimated parameters X:
%        Y = f(X5,LMi,LMf);
        
  %%%% After getting the estimated parameters, calculate the phasors and
  %%%% errors
  %%%%% Intermediate phasors calculation
        T = NSamples*dt;
%        tau_est2 = 0.25;
        tau_est5 = tau_pp + tau/NSamples + 2*dt/T;  %estimated tau/T (we can simulate errors in tau estimation)
%        tau_est7 = 0.75;
        if ha ~= 0
            %phase step
            Xe5_r = xr(1)/sqrt(2);   %reference magnitude (actual value xr)
            Xe5 = X5(1)/sqrt(2);      %estimated magnitude
            Phe_r = xr(3)*tau_pp + (xr(3) + xr(4)*(pi/180))*(1 - tau_pp); %reference phase (actual value xr)
            Phe5 = X5(3)*tau_est5 + (X5(3)+X5(4)*(pi/180))*(1 - tau_est5);  %estimated phase
        else
            %mag step
%            Xe2_r = (xr(1)*tau_est2 + xr(1)*(1+xr(2))*(1 - tau_est2))/sqrt(2);  %reference magnitude
            Xe5_r = (xr(1)*tau_est5 + xr(1)*(1+xr(2))*(1 - tau_est5))/sqrt(2);
%            Xe7_r = (xr(1)*tau_est7 + xr(1)*(1+xr(2))*(1 - tau_est7))/sqrt(2);
                        
%            Xe2 = (X2(1)*tau_est2 + X2(1)*(1+X2(2))*(1-tau_est2))/sqrt(2);
            Xe5 = (X5(1)*tau_est5 + X5(1)*(1+X5(2))*(1-tau_est5))/sqrt(2); %estimated magnitude
%            Xe7 = (X7(1)*tau_est7 + X7(1)*(1+X7(2))*(1-tau_est7))/sqrt(2);
            
            Phe_r = xr(4); %reference phase [rad]
            
%            Phe2 = X2(4);% - 0.25*(X2(3)-F0*2*pi); correction?
            Phe5 = X5(4);    %estimated phase
%            Phe7 = X7(4);% + 0.25*(X7(3)-F0*2*pi);
        end


%%%%% Adjacent window calculation
    SignalParams(1) = F1; SignalParams(2:8) = 0;
    DelayCorr = 0; %-3500.0; %   %[in nanosecond]
    MagCorr = 1;
    [Synx_bw,Freq_bw,ROCOF_bw] = SteadyStateFit ( ...
	SignalParams, ...
	DelayCorr, ...
	MagCorr, ...
	F0, ...
	AnalysisCycles, ...
	SampleRate, ...
	Signal(NSamples/2-499:NSamples/2) ...
);
        
[Synx_aw,Freq_aw,ROCOF_aw] = SteadyStateFit ( ...
	SignalParams, ...
	DelayCorr, ...
	MagCorr, ...
	F0, ...
	AnalysisCycles, ...
	SampleRate, ...
	Signal(NSamples/2+1:NSamples/2+500) ...
);


        %TVE
        TVE = @(Ve_r,Ve) 100*sqrt(((real(Ve_r) - real(Ve))^2 + (imag(Ve_r) - imag(Ve))^2 )/(real(Ve_r)^2 + imag(Ve_r)^2)) ;

        %rectangular Reference phasors
        corr = (F1-F0)*2*pi*(250/5000)/2;
        Ve_ref(1:10) = (xr(1)/sqrt(2))*exp(1i*Phe_r);
        Ve_ref(4) = Xe5_r.*exp(1i*Phe_r);        
        if ha ~= 0
            Ve_ref(6:10) = xr(1)*exp(1i*(Phe_r + xr(4)));
        else
            Ve_ref(6:10) = ((xr(1)*(1+xr(2)))/sqrt(2))*exp(1i*Phe_r);
        end
        
        %Ve2 = Xe2*exp(1i*(Phe2 + corr));
        Ve5 = Xe5*exp(1i*Phe5);
        %Ve7 = Xe7*exp(1i*(Phe7 - corr));
        Ve(1:3) = Synx_bw;
        Ve(4) = [Ve5];
        Ve(5:10) = Synx_aw;
        Freq_LM(1:3) = Freq_bw;
        if ha ~= 0
            Freq_LM(4) = [X5(2)/(2*pi)];
        else 
            Freq_LM(4) = [X5(3)/(2*pi)];
        end
        Freq_LM(5:10) = Freq_aw;
        Freq_SS(1:3) = Freq_bw;
        Freq_SS(4) = (Freq_bw+Freq_aw)/2;
        Freq_SS(5:10) = Freq_aw;
        
        Synx_ref = Ve;
        Synx_ref(3:4) = Synx_bw;
        Synx_ref(5) = Synx_aw;
        
        %intermediate phasors from SS
        Ve_synx(1:10) = Synx_bw;
        %Ve_synx(3) = (0.75*abs(Synx_bw)+0.25*abs(Synx_aw))*exp(1i*(0.75*angle(Synx_bw) + 0.25*angle(Synx_aw)));
        Ve_synx(4) = (0.5*abs(Synx_bw)+0.5*abs(Synx_aw))*exp(1i*(0.5*angle(Synx_bw) + 0.5*angle(Synx_aw)));
        %Ve_synx(5) = (0.25*abs(Synx_bw)+0.75*abs(Synx_aw))*exp(1i*(0.25*angle(Synx_bw) + 0.75*angle(Synx_aw)));
        Ve_synx(6:10) = Synx_aw;

        %Freq_ref
        Freq_ref(1:10) = F1;
        Freq_ref(4) = F1;

        %PMU phasor obtained by fft
%         ini = 501; m = 250; wsize = 999;

        for pm = 1:10
%             PMU_signal = Signal(ini + pm*m:ini + pm*m + wsize);
%             %plot(PMU_signal)
%             %PMU_signal = Signal(ini+m:1750);
%             [V_pmu(pm), pmu_freq(pm)] = PMU_FFT(PMU_signal,SampleRate);
%             %PMU/DFT comparison
%             TVE_aw(pm) = TVE(Synx_ref(pm),V_pmu(pm));
%             TVE_ir(pm) = TVE(Ve(pm),V_pmu(pm));
%             
%             %PMU/DFT with intermediate phasors from SS
%             TVE_irs(pm) = TVE(Ve_synx(pm),V_pmu(pm));
            
            %comparison of estimators
            TVE_refLM(pm) = TVE(Ve_ref(pm),Ve(pm));
            TVE_refSS(pm) = TVE(Ve_ref(pm),Ve_synx(pm));
            
            %comparison of frequency
            FE_refSS(pm) = abs(Freq_SS(pm) - Freq_ref(pm));
            FE_irs(pm) = abs(Freq_LM(pm) - Freq_ref(pm));
           
        end
%         tt = (1:10)*m/SampleRate;
%         stdlimit = ones(1,10);
%         subplot(2,1,1)
%         plot(tt,stdlimit,'k',tt,TVE_aw,'k^--',tt,TVE_ir,'ro:'); ylabel('TVE[%]');xlabel('Time [s]');
%         ylim([-0.5 5.5])
%         title('TVE of a DFT PMU');
%         legend('Standard limit','Reference from adjacent windows','Reference from intermediate phasors')
%         subplot(2,1,2)
%         plot(tt,TVE_ir,'ro:',tt,TVE_irs,'bx--'); ylabel('TVE[%]');xlabel('Time [s]');
%         legend('Reference from LM','Reference from adjacent windows')
%         title('TVE of a DFT PMU with intermediate reference phasors');
 
%        figure
%        plot(1:10,TVE_refLM,'ro:',1:10,TVE_refSS,'bx--'); ylabel('TVE[%]');xlabel('Time [s]');
        %ylim([1e-6 1e-1])
        %legend('Reference from LM','Reference from SS adjacent windows')
        %title('TVE of estimators with intermediate reference phasors');
        
%        figure
%        plot(1:10,FE_irs,'ro:',1:10,FE_refSS,'bx--'); ylabel('FE[Hz]');xlabel('Time [s]');
%        legend('F_{LM} - F1','F_{SS} - F1')
%        title('FE of estimators with intermediate reference phasors');

        Ve_s_ang = angle(Ve_synx)*180/pi;
        Ve_LM_ang = angle(Ve)*180/pi;
        Ve_ref_ang = angle(Ve_ref)*180/pi;
        
        Ph_errorSS = Ve_s_ang(3:5) - Ve_ref_ang(3:5);
        Ph_errorLM = Ve_LM_ang(3:5) - Ve_ref_ang(3:5);

%         figure
%         plot(tt(3:5),Ph_errorSS,'o:',tt(3:5),Ph_errorLM ,'x--')
%         title('Phase errors');legend('SS','LM')
        
        mag_SS = abs(Ve_synx);
        mag_LM = abs(Ve);
        mag_ref = abs(Ve_ref);
        Mag_errorSS = mag_SS - mag_ref;
        Mag_errorLM = mag_LM - mag_ref;
        
        FE_LM = Freq_LM(4) - F1;
        FE_SS = Freq_SS(4) - F1;
        
        TVE_LM = TVE_refLM(4);
        TVE_SS = TVE_refSS(4);
       
        