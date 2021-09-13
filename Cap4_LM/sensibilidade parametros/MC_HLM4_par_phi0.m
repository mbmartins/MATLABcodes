function [filename, std_R2, std_w, std_Xmag, std_Xphi] = MC_HLM4_par_phi0(MCruns, par_var, SNR, Fs, hm, ha, phi0,run)
% só há variação em phi0, os outros parms ficam constantes
% Estimation of Sincrophasors with discontinuities
% using MATLAB optimization toolbox 

%signal generation
F0 = 60; %Hz
F1 = 60; %Hz
%Fs = 5000; %Hz
dt = 1/Fs;
AnalysisCycles = 6;
NSamples = floor(AnalysisCycles*Fs/F0);
n = -NSamples/2:(NSamples/2-1); %discrete time vector

Xm = 1; %70*sqrt(2) =~ 100;
%phi0 = -120; %phase in degrees
Ph = phi0*pi/180;% Phase in radians
%ha = -10;   % IEEE Std phase (angle) step index: 10 degrees
%hm = 0;   % magnitude step index: 0.1 
Wf = 2*pi*F1;  % fundamental frequency
%SNR = 60; %dB SNR = 20 log_10 Asinal/Aruido => Aruido = Asinal/10^(SNR/20)
Aruido = Xm/10^(SNR/20);

if run == true
    for ti = 1:9

        tau_pp = 0.1*ti; % relative time of step in percent of total time 
        tau_0 = (tau_pp - 0.5)*NSamples; %discrete time displacement
        n = -NSamples/2:(NSamples/2-1); %discrete time vector
        n = n - tau_0;
        t = n*dt; %time vector

        %%%% Signal parameters estimation

        u = zeros(length(Xm),length(t));
        tau = 0;
        u(length(Xm),t >= tau) = u(length(Xm),t >= tau) + 1;

        %modelo 1
        % Phase step
        % Model: f(x) = x1*cos(w*t + phi + x2*u(t - tau)) 
        if ha ~= 0
            %f = @(x) x(1)*cos(x(2)*t + x(3) + x(4)*(pi/180)*u);
            xnom = [Xm 2*pi*F1 Ph ha];
            f = @(x) xnom(1)*cos(xnom(2)*t + x(3) + 2*pi + xnom(4)*(pi/180)*u);
        else
            %mag
            %f = @(x) x(1)*(1+x(4)*u).*cos(x(2)*t + x(3));
            %xnom1 = [Vm 2*pi*F1 Ph KxS];
            xnom = [Xm hm 2*pi*F1 Ph];
            f = @(x) xnom(1)*(1+xnom(2)*u).*cos(xnom(3)*t + x(4));

        end

        %err = @(x) (Signal - f(x)).^2;
        err = @(x) (Signal - f(x));

        % Nonlinear fit
        % Monte Carlo analysis
        %MCruns = 1000;
        x0 = xnom;  %x0 is fixed - first guess are the nominal values

        for k = 1:MCruns
            %first guess

            %uncertainties of parameters in signal generation
%             if ha ~= 0
%             % phase         X1  w    ph   x3 (KaS)
%                 par_var = [1   0.05 1 1]; % parameter variation in percent related to nominal
%             else
%             % mag          x1  x2(KxS)  wf    ph  
%                 par_var = [1   1     0.05  1]; % parameter variation in percent related to nominal            
%             end

            rng('shuffle');
            rn = (rand(1,length(xnom))-0.5);
            xr = xnom.*(1+2*(par_var/100).*rn);
            freq_rand = xr(3)/(2*pi);

            %uncertainties of tau estimation
            utau = 0;  %number of dts 
            u = zeros(length(Xm),length(t));
            tau = dt*randi([-utau utau],1);
            u(length(Xm),t >= tau) = u(length(Xm),t >= tau) + 1;

            %regenerates the signal with the uncertainties
                i = 1;
                Xm = xr(1); 
                if ha ~= 0
                    Wf(i) = xr(2); 
                    Ph(i) = xr(3);
                    KxS_(i) = 0;
                    KaS_(i) = xr(4);
                else
                    Ph(i) = xr(4); 
                    Wf(i) = xr(3);
                    KxS_(i) = xr(2);
                    KaS_(i) = 0;
                end

                Ain = zeros(length(Xm),length(t));
                % Amplitude Step: applied after time passes 0
                Ain(i,:) = Xm(i);
                Ain(i,t >= 0) = Ain(i,t >= 0) * (1 + KxS_(i));
                %Phase step
                Theta(i,:) = (Wf(i)*t) ...                         % Fundamental
                                 + Ph(i);               % phase shift
                Theta(i,t >= 0) = Theta(i,t >= 0) + (KaS_(i) * pi/180);
                cSignal = (Ain.*exp(-1i.*Theta));

                rSignal = real(cSignal);
                var_noise = ((std(rSignal))/(10^(SNR/20)))^2;
                std_noise = (std(rSignal))/(10^(SNR/20));
                noise = std_noise*randn(1,length(rSignal));
                Signal = rSignal + noise;
                SNR_hist = snr(Signal,noise);
                %err = @(x) (Signal - f(x)).^2;
                err = @(x) (Signal - f(x));
            y0 = f(x0);
            x = xnom;
            %estimated signal
    %          y = f(xnom);
    %          plot(t,Signal,'.b',t,y,'r')

            % figure
            % plot(t,Signal - y0,'r')
            % title('Point by point estimation error: Signal - y')

            % Levenberg-marquardt from optimization toolbox
            tol = 1e-7;
            OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','OptimalityTolerance',tol);
            OPTIONS.StepTolerance = 1e-12;
            OPTIONS.Display = 'off';
            [X,RESNORM(k),RESIDUAL,exitflag,output] = lsqnonlin(err,x0,[],[],OPTIONS);
            Y = f(X);
            %RESNORM = R2(k) = sum(RESIDUAL.^2);
                      w(k) = X(3);
        %Intermediate phasor
            T = NSamples*dt;
            tau_est = tau_pp + tau/NSamples;            
        if ha ~= 0 % salto fase
              w(k) = X(2);
            Xmag(k) = X(1);
            Xphi(k) = X(3)*tau_est + (X(3)+X(4)*(pi/180))*(1 - tau_est);
        else
               w(k) = X(3);
            Xmag(k) = X(1)*tau_est + X(1)*(1+X(2))*(1-tau_est);
            Xphi(k) = X(4);
        end           
            fprintf("ti = "+ ti + "; k = " + k + "; R2 = "+ RESNORM(k) + "\n");
        end
        std_R2(ti) = std(RESNORM);
        std_w(ti) = std(w);    
        std_Xmag(ti) = std(Xmag);
        std_Xphi(ti) = std(Xphi);
    end    
end


filename = "Cap4_LM\resultadosMC\MC_LM_RESIDUOS_";
par_str = "PAR"+par_var(1)+par_var(2)+par_var(3)+par_var(4);
filename = filename+par_str+"MCruns"+MCruns+"Fs"+Fs+"SNR"+SNR+"ha"+ha+"hm"+(hm*100)+"phi0"+phi0
save(filename)

% if KaS ~= 0
%     plot(errors(:,3),errors(:,4),'.')
%     title('Phase Step: Correlation \phi vs X3')
% else
%     plot(errors(:,1),errors(:,2),'.')
%     title('Mag Step: Correlation \epsilon _{x_1} vs \epsilon _{x_2}')
% end
