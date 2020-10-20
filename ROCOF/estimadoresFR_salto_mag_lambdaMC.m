%testes para estimação de fase, frequência e ROCOF

%variacoes interessantes:
% f1 = 59, 61, etc
% dois saltos
% encontrar um lambda otimo para o PATV

clear all; close all; clc
SNR = 60;
%fixed parameters
F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency
KaS = 0.0; %[degrees]
KxS = 0.1; % [relative magnitude step]
KfS = 0.; %[Hz] %size of frequency step
Fs =4800;
NCycles = 6;
T = NCycles/F0;
phi_0 = 0; %angle phi_0 [degrees]
tau1 = 0.5;  % in proportion of T
tau2 = 1.; % in proportion of T; set tau2 = 1 if you dont want two steps
NSamples = floor(NCycles*Fs/F0);
tau_n1 = floor(tau1*NSamples);
tau_n2 = floor(tau2*NSamples);
nbits = 16;

% vetor de valores de lambda
lambda_step = .05;
n_lambdas = 100; la_ini = .001;
lambda_n = la_ini+(0:lambda_step:(n_lambdas-1)*lambda_step);
% --- grade fina, valores pequenos de lambda
% lambda_step = .01;
% n_lambdas = 100; la_ini = 0.01;
% lambda_n = la_ini+(0:lambda_step:(n_lambdas-1)*lambda_step);

%teste de lambda unico
%lambda_n = 5.48;n_lambdas = 1;

MCiter = 300;
phi_n = 360*rand(1,MCiter); %distribuicao de phi_0
%tau_vec = 0.5*ones(1,MCiter);
%tau_vec = 0.1:0.1:0.5;
tau_vec = 0.1 + (0.8)*rand(1,MCiter); %distribuicao de tau

for L=1:n_lambdas
    L
    % MC loop
    for k=1:MCiter
        Signal = SigGEN2(F0,F1,Fs,phi_n(k),NCycles,tau_vec(k),tau2,SNR,KaS, KxS,KfS,nbits);
        z=hilbert(Signal');
        Psi_i = unwrap(angle(z)); %[rad]
        %Psi_i = phase(z);
        f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
        az = abs(z);

        tau = tau_vec(k)*T;
        Fref = (tau*F1 + (T-tau)*(F1 + KfS))/T;
        tau_n1 = floor(tau*NSamples/T);

        %estimador EFx
        [f1_est4(k),f2_est4(k),F_est4(k),fu,ri] = EF4(Psi_i,az,Fs,tau_n1,lambda_n(L));
        [f1_est5(k),f2_est5(k),F_est5(k),fu,ri] = EF5(f_i,az,tau_n1,lambda_n(L));
        [f1_est6(k),f2_est6(k),F_est6(k),fu,ri] = EF6(f_i,az,tau_n1,lambda_n(L));
        
        % ---- DEBUG ----
%         F_mean(k) = mean(F_est); 
%         fig1 = figure(1)
%          plot(F_mean,'x-');
%          pause(0.5)
    end
        FE4(L,:) = F_est4 - Fref;
        FE14(L,:) = f1_est4 - F1;
        kfE4(L,:) = f2_est4 - f1_est4 - KfS;    

        FE5(L,:) = F_est5 - Fref;
        FE15(L,:) = f1_est5 - F1;
        kfE5(L,:) = f2_est5 - f1_est5 - KfS;
        
        FE6(L,:) = F_est6 - Fref;
        FE16(L,:) = f1_est6 - F1;
        kfE6(L,:) = f2_est6 - f1_est6 - KfS;
       
        
end

save('salto_mag_lambdaMC')
% se desejar figuras rodar
run('analise_lambda.m')


% FEmax4 = max(abs(FE4'));
% FEmax5 = max(abs(FE5'));
% FEmax6 = max(abs(FE6'));
% 
% FE1max6 = max(abs(FE16'));
% % FE2max(L,tau_i) = max(abs(FE2));
% % kfEmax(L,tau_i) = max(abs(kfE));
% 
% FEmean4 = mean(FE4');
% FEmean5 = mean(FE5');
% FEmean6 = mean(FE6');
% FE1mean6 = mean(FE16');
% 
% 
% FE_std4 = std(FE4');
% FE_std5 = std(FE5');
% FE_std6 = std(FE6');
% 
% %----- Figuras para FE x lambda
% figs(1) = figure(1); hold off
% semilogy(lambda_n,abs(FEmean4),'b.-'); hold on;
% semilogy(lambda_n,abs(FEmean4)+FE_std4,'k--');
% semilogy(lambda_n,(FEmax4),'k:');
% grid on;
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$|FE|$ [Hz]')
% legend('|\mu(FE)|','|\mu(FE)| + \sigma(FE)', '|FE_{max}|')
% title('EF4')
% 
% figs(2) = figure(2);hold off
% semilogy(lambda_n,abs(FEmean5),'b.-'); hold on;
% semilogy(lambda_n,abs(FEmean5)+FE_std5,'k--');
% semilogy(lambda_n,(FEmax5),'k:');
% grid on;
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$|FE|$ [Hz]')
% legend('|\mu(FE)|','|\mu(FE)| + \sigma(FE)', '|FE_{max}|')
% title('EF5')
% 
% figs(3) = figure(3); hold off
% semilogy(lambda_n,abs(FEmean6),'b.-'); hold on;
% semilogy(lambda_n,abs(FEmean6)+FE_std6,'k--');
% semilogy(lambda_n,(FEmax6),'k:');
% grid on;
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$')
% ylabel('$|FE|$ [Hz]')
% legend('|\mu(FE)|','|\mu(FE)| + \sigma(FE)', '|FE_{max}|')
% title('EF6')
% 
% savefig(figs,'salto_fase_lambda')
% save('salto_fase_lambda_MC')
% 
% % ---- Figuras para FE1 x lambda
% % figure
% % semilogy(lambda_n,abs(FE1mean6),'o-')
% % xlabel('x','Interpreter','latex');
% % ylabel('y','Interpreter','latex');
% % xlabel('$\lambda$')
% % ylabel('$FE1_{\mu}$ [Hz]')
% % 
% % 
% % figure
% % plot(lambda_n,FE2max,'o-')
% % xlabel('x','Interpreter','latex');
% % ylabel('y','Interpreter','latex');
% % xlabel('$\lambda$')
% % ylabel('$FE2_{max}$ [Hz]')
% % legend('\tau = 0.1T','\tau = 0.2T', '\tau = 0.3T','\tau = 0.4T','\tau = 0.5T')