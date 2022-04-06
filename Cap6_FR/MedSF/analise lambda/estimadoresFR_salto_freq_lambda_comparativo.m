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
phi_0 = 100; %angle phi_0 [degrees]
tau1 = 0.5;  % in proportion of T
tau2 = 1.; % in proportion of T; set tau2 = 1 if you dont want two steps
N = floor(NCycles*Fs/F0);
tau_n1 = floor(tau1*N); 
% tau_n2 = floor(tau2*N);
% nbits = 16;

lambda_n = [0.5 3 30]

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
    
        figure(1)
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
        [f1_est(3),f2_est(3),F_est(3),fu,ri] = MedSF_PATV(f_i,az,tau_n1,lambda_n(3));
        hf_est = f2_est(3) - f1_est(3);
        hfE(3) = hf_est - hf;
        FE(3) = F_est(3) - Fref;
        %plot(n(brmask),f_i(brmask)); hold on;
        plot(brmask,fu,'g')
        lg = legend('$f_i[n] [Hz]$','$s[n]+c_f[n]$, $\lambda = 0,5$','$s[n]+c_f[n]$, $\lambda = 3$','$s[n]+c_f[n], \lambda = 30$')
        lg.Interpreter = 'latex';lg.FontSize = 13;
        grid on;
        ylb = ylabel('$f_i[n]$'); ylb.Interpreter = 'latex';ylb.FontSize = 13;
        xlb = xlabel('$n$'); xlb.Interpreter = 'latex';xlb.FontSize = 13;
        
        %saveas(gcf,'efeito_lambda_decomposicao_fi.png')
        
        figure(2)
        pl = line([q,tau_n1],[F1,F1])
        pl1 = line([q,tau_n1],[f1_est(1),f1_est(1)]);pl1.Color = 'red'; pl1.LineStyle = '--'; % f1 estimado
        pl2 = line([q,tau_n1],[f1_est(2),f1_est(2)]);pl2.Color = 'black'; pl2.LineStyle = '--'; % f1 estimado
        pl3 = line([q,tau_n1],[f1_est(3),f1_est(3)]);pl3.Color = 'green'; pl3.LineStyle = '--'; % f1 estimado
        F2 = F1 + hf;
        pl = line([tau_n1,N-q],[F2,F2])
        pl1 = line([tau_n1,N-q],[f2_est(1),f2_est(1)]);pl1.Color = 'red'; pl1.LineStyle = '--'; % f2 estimado
        pl2 = line([tau_n1,N-q],[f2_est(2),f2_est(2)]);pl2.Color = 'black'; pl2.LineStyle = '--'; % f2 estimado
        pl3 = line([tau_n1,N-q],[f2_est(3),f2_est(3)]);pl3.Color = 'green'; pl3.LineStyle = '--'; % f2 estimado
        lg = legend('$f[n]$', '$\hat{f}_1, \hat{f}_2, \lambda = 0,5$','$\hat{f}_1, \hat{f}_2,\lambda = 3$','$\hat{f}_1, \hat{f}_2,\lambda = 30$')

        lg.Interpreter = 'latex';lg.FontSize = 13;
        grid on;
        ylb = ylabel('$f_i[n]$'); ylb.Interpreter = 'latex';ylb.FontSize = 13;
        xlb = xlabel('$n$'); xlb.Interpreter = 'latex';xlb.FontSize = 13;
                saveas(gcf,'efeito_lambda_decomposicao_f1_f2.png')