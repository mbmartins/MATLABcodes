%testes para estimação de fase, frequência e ROCOF

% encontrar um lambda otimo para o PATV

clear all; close all; clc
SNR = 100;
%fixed parameters
F0 = 60.0; %nominal frequency
F1 = 58; %fundamental frequency
KaS = 0.0; %[degrees]
KxS = 0.0; % [relative magnitude step]
hf = -1.; %[Hz] %size of frequency step
Fs =4800;
NCycles = 6;
T = NCycles/F0;
phi_0 = 0; %angle phi_0 [degrees]
tau1 = 0.5;  % in proportion of T
tau2 = 1.; % in proportion of T; set tau2 = 1 if you dont want two steps
N = floor(NCycles*Fs/F0);
tau_n1 = floor(tau1*N); 
% tau_n2 = floor(tau2*N);
% nbits = 16;

lambda_n = [0.5 3]

n=0:N-1; n=n(:); %base de tempo

L=1;

tau = tau1*T;
        
Fref = (tau*F1 + (T-tau)*(F1 + hf))/T;
tau_n1 = round(tau*N/T);
        
% --- gerador sinal CA com salto frequencia
        w1 = 2*pi*F1/Fs; w2 = 2*pi*hf/Fs; Xm = 1;
        phi_0_rad=phi_0*pi/180; % phi_0 sorteado
        phi_0_rad = phi_0_rad + 2*pi*(F0 - F1)*tau_n1/Fs; % fator de correção para F1 fora da nominal
        PHI=w1*n+w2.*(n-tau_n1).*(n>=tau_n1)+phi_0_rad; % fase instantanea
        x=Xm.*cos(PHI);  % sinal x[n] 
        vx=var(x); vruido=vx./(10^(SNR/10)); 
        xn=x+sqrt(vruido)*randn(size(n));  % sinal AC ruidoso]
          
        z=hilbert(xn);
        zx = hilbert(x);
        %Psi_i = unwrap(angle(z)); %[rad]
        Psi_i = phase(z);
        Psi_i_x = phase(zx);
        f_ix = gradient(Psi_i_x)*Fs/(2*pi);
        f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
        az = abs(z);
        %estimadores
                
        
        
        n = 1:N;
    br = 0.05;%80/480; % fraction of samples to be ignored 
    q = floor(br*N); % number of samples to be ignored
    brmask = [q+1:N-q];
        [f1_est(1),f2_est(1),F_est(1),fu,ri] = MedSF_PATV(f_i,az,tau_n1,lambda_n(1));
        FE(1) = F_est(1) - Fref;
        hf_est = f2_est(1) - f1_est(1);
        hfE(1) = hf_est - hf; 
        plot(n(brmask),f_i(brmask),'b:'); hold on;
        plot(brmask,fu,'r')
        [f1_est(2),f2_est(2),F_est(2),fu,ri] = MedSF_PATV(f_i,az,tau_n1,lambda_n(2));
        hf_est = f2_est(2) - f1_est(2);
        hfE(2) = hf_est - hf;
        FE(2) = F_est(2) - Fref;
        %plot(n(brmask),f_i(brmask)); hold on;
        plot(brmask,fu,'k')
        grid on;
        plot(brmask,f_ix(brmask),'b')
        lg = legend('$f_i[n]$','$\hat{f}_i[n]$','$s[n]+c_f[n]$ $\lambda = 0,5$','$s[n]+c_f[n]$ $\lambda = 3$')
        lg.Interpreter = 'latex';
        lg.FontSize = 13;
        xlb = xlabel('n'); 
        ylb = ylabel('$\hat{f}_i[n]$ [Hz]');
        xlb.Interpreter = 'latex';
        xlb.FontSize = 13;
        ylb.Interpreter = 'latex';
        ylb.FontSize = 13;

        [f1_est(3),f2_est(3),F_est(3),fu,~] = MedSF(f_i, tau_n1);
        hf_est = f2_est(3) - f1_est(3);
        hfE(3) = hf_est - hf;
        FE(3) = F_est(3) - Fref;