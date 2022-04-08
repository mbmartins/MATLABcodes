%gera grafico do teste de duracao da janela
clear all; close all; clc;
F0 = 60;
F1_vec = [59.5 60.0 60.5 61];
Fs = 4800;
hf = -2;
phi_n = 0;
tau_n = 240;
NA = 480;
n = 1:NA;
q = round(0.05*NA);

%set(groot,'defaultAxesColorOrder',[0 0 1;0 0 0.75;0 0 0.25;0 0 0]);

for m=1:length(F1_vec)
        F1 = F1_vec(m);
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
        figure(1); colormap winter;
        pl = plot(f_i(q:NA-q)); hold on; %plot(tau_nB,Xm,'x')
        %axis([1 NA -1 1])
        xlabel('$n$','Interpreter','latex')
        ylabel('$f_i[n]$','Interpreter','latex')
        %legend('$f_1='+F1+'$')
        grid on
end
lg = legend('$f_1=59,5$','$f_1=60,0$','$f_1=60,5$','$f_1=61$');
lg.Interpreter='latex';
set(groot,'defaultAxesLineStyleOrder','remove')
set(groot,'defaultAxesColorOrder','remove')