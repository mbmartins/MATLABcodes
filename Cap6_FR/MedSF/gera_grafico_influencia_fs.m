%gera grafico do teste de duracao da janela
clear all; close all; clc;
F0 = 60;
F1_vec = [60];
Fs_vec = [2400 9600 192000];
hf_vec = -1;
phi_n = 0;
T = 0.1;
NA = T*Fs_vec;



%set(groot,'defaultAxesColorOrder',[0 0 1;0 0 0.75;0 0 0.25;0 0 0]);

for m=1:length(Fs_vec)
        n = 1:NA(m);
        q = round(0.05*NA(m));
        tau_n = NA(m)/2;
        F1 = F1_vec;
        hf = hf_vec;
        Fs = Fs_vec(m);
        t = 0:1/Fs:T;
        fsf1 = Fs/F1;
    %         % --- gerador 
        w1 = 2*pi*F1/Fs; w2 = 2*pi*hf/Fs; Xm = 1;
        phi_0_rad=phi_n*pi/180; % phi_0 sorteado
        phi_0_rad = phi_0_rad + 2*pi*(F0 - F1)*tau_n/Fs; % fator de correção para F1 fora da nominal
        PHI=w1*n+w2.*(n-tau_n).*(n>=tau_n)+phi_0_rad; % fase instantanea
        x=Xm.*cos(PHI);  % sinal x[n] 
        %subplot(3,1,m)
        z = hilbert(x);
        Psi_i = phase(z);
        f_i=gradient(Psi_i)*Fs/(2*pi);
        figure(1); %colormap winter;
        %subplot(211)
        pl = plot(t(q:NA(m)-q),f_i(q:NA(m)-q)); hold on; %plot(tau_nB,Xm,'x')
        pl.Marker = '.';
        %axis([1 NA -1 1])
        xlabel('$t$ [s]','Interpreter','latex')
        ylabel('$f_i[n]$','Interpreter','latex')
        %legend('$f_1='+F1+'$')
        f1_est(m) = median(f_i(q:tau_n-1));
        f2_est(m) = median(f_i(tau_n+1:NA(m)-q));

        grid on
%         figure(2); hg = histogram(f_i(q:tau_n-1)); hold on;
%         hg.Normalization = 'probability';
%         hg.BinLimits = [60 61];
%         hg.NumBins = 13;
end
lg = legend('$f_s=2400$','$f_s=9600$','$f_s=19200$');
figure(1)
lg.Interpreter='latex';
%subplot(212)


% F2 = F1 + hf;
%         FE1 = f1_est - F1;
%         FE2 = f2_est - F2;
%         line([0,0.05],[FE1(1),FE1(1)],'Color','blue');
%         line([0.05,0.1],[FE2(1),FE2(1)],'Color','blue');
%                 line([0,0.05],[FE1(m),FE1(m)],'Color','red');
%         line([0.05,0.1],[FE2(m),FE2(m)],'Color','red');
%         lg = legend('$f_s=2400$','$f_s=2400$','$f_s=19200$','$f_s=19200$')
% grid on