%testes para estima��o de fase, frequ�ncia e ROCOF

% encontrar um lambda otimo para o PATV

clear all; close all; clc
SNR = 60;
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

% vetor de valores de lambda
 lambda_step = .02;
 n_lambdas = 100; 
 la_ini = .01;
 lambda_n = la_ini+(0:lambda_step:(n_lambdas-1)*lambda_step);

% --- grade fina, valores pequenos de lambda
% lambda_step = .01;
% n_lambdas = 100; la_ini = 0.01;
% lambda_n = la_ini+(0:lambda_step:(n_lambdas-1)*lambda_step);

%teste de lambda unico
%lambda_n = 0.5;n_lambdas = 1;

MCiter = 300;
%phi_0v = 180*rand(1,MCiter); %distribuicao de phi_0
%tau_vec = 0.1 + (0.8)*rand(1,MCiter); %distribuicao de tau
%ert_vec = randi([-8 8],1,MCiter); % erros de tau_n de -8 a 8

% valores fixos
 phi_0v = 0*ones(1,MCiter);
 tau_vec = 0.5*ones(1,MCiter);
ert_vec = zeros(1,MCiter);



n=0:N-1; n=n(:); %base de tempo

for L=1:n_lambdas
    %L
    % MC loop
    for k=1:MCiter
        %sorteio de phi0 e tau
        phi_0 = phi_0v(k);
        tau1 = tau_vec(k);
        erro_tau = ert_vec(k);

%         %fixos para debug
%           phi0 = 0; % em graus
%           tau1 = 0.1;
        
        tau = tau1*T;
        
%        Signal = SigGEN2(F0,F1,Fs,phi_0,NCycles,tau1,tau2,SNR,KaS, KxS,hf,nbits);

        Fref(k) = (tau*F1 + (T-tau)*(F1 + hf))/T;
        %tau_n1 = floor(tau*NSamples/T);
        tau_n1 = round(tau*N/T);
        
%         % --- gerador alternativo
        w1 = 2*pi*F1/Fs; w2 = 2*pi*hf/Fs; Xm = 1;
        phi_0_rad=phi_0*pi/180; % phi_0 sorteado
        phi_0_rad = phi_0_rad + 2*pi*(F0 - F1)*tau_n1/Fs; % fator de corre��o para F1 fora da nominal
        PHI=w1*n+w2.*(n-tau_n1).*(n>=tau_n1)+phi_0_rad; % fase instantanea
        x=Xm.*cos(PHI);  % sinal x[n] 
        vx=var(x); vruido=vx./(10^(SNR/10)); 
        xn=x+sqrt(vruido)*randn(size(n));  % sinal AC ruidoso]
        
        %plot(xn);hold on; plot(Signal)
    
        z=hilbert(xn);
        %Psi_i = unwrap(angle(z)); %[rad]
        Psi_i = phase(z);
        f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
        az = abs(z);
        %estimador
        tau_n1e = tau_n1 + erro_tau;
        [f1_est(k),f2_est(k),F_est(k),fu,ri] = MedSF_PATV(f_i,az,tau_n1e,lambda_n(L));
        %sprintf('%d',k)
    end

        FE(L,:) = F_est - Fref;
        FE1(L,:) = f1_est - F1;
        FE2(L,:) = f2_est - F1 - hf;
        hfE(L,:) = f2_est - f1_est - hf;
        hFEmean = abs(mean(hfE'));
        str = "lambda=" + sprintf('%d',lambda_n(L)) + "; HFEmean =" + sprintf('%f',hFEmean(L))
        figure(1)
        semilogy(lambda_n(1:length(hFEmean)),abs(hFEmean),'k.-');
end

%save('salto_freq_lambdaMC_MedSF300_v2')

FE_std = std(FE');
FEmean = abs(mean(FE'))

FE1_std = std(FE1');
FE1mean = abs(mean(FE1'));

FE2_std = std(FE2');
FE2mean = abs(mean(FE2'));

HFE_std = std(hfE');
HFEmean = abs(mean(hfE'));

close all;
% %----- Figuras para FE x lambda
% 
% figs(2) = figure(2);hold off
% 
% %subplot(131)
% semilogy(lambda_n,abs(FEmean),'b.-');hold on; 
% semilogy(lambda_n,abs(FEmean) + FE_std,'b-.');
% semilogy(lambda_n,abs(FE1mean),'r.-'); hold on;
% semilogy(lambda_n,abs(FE1mean)+FE1_std,'r-.');
% semilogy(lambda_n,abs(FE2mean),'g.-'); hold on;
% semilogy(lambda_n,abs(FE2mean)+FE2_std,'g-.');
% semilogy(lambda_n,abs(HFEmean),'k.-'); hold on;
% semilogy(lambda_n,abs(HFEmean)+HFE_std,'k-.');

% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% xlabel('$\lambda$','Fontsize',16)
% ylabel('$|FE|,|FE_1|,|h_fE|$ [Hz]','Fontsize',15)
% legend('|\mu(FE)|','|\mu(FE)| + \sigma(FE)',...
%        '|\mu(FE_1)|','|\mu(FE_1)| + \sigma(FE_1)',...
%        '|\mu(FE_2)|','|\mu(FE_2)| + \sigma(FE_2)',...
%        '|\mu(h_fE)|','|\mu(h_fE)| + \sigma(h_fE)')
% legend('Location','eastoutside')
% legend('Fontsize',12)
% %xlabel('Fontsize',14)
% %ylabel('Fontsize',14)
% legend('Orientation','vertical')
% grid on
% %ylim([10e-5 1]);
% %xlim([0.1 7])

save('lambda_erro_tau.mat')
run('analise_lambda_corrigido')