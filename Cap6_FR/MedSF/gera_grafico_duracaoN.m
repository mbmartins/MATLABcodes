%gera grafico do teste de duracao da janela
clear all; close all; clc;
F0 = 60;
F1 = 60;
Fs = 4800;
hf = -1;
phi_n = 0;
tau_nA = 240;
NA = 480;
fsf1 = Fs/F1;

for m = 1:3
    tau_nB(m) = tau_nA + (m-2)*fsf1;
    NB(m) = round(tau_nB(m)*(NA/tau_nA));
end

for m=1:3
%         % --- gerador 
    n = 1:NB(m);
w1 = 2*pi*F1/Fs; w2 = 2*pi*hf/Fs; Xm = 1;
        phi_0_rad=phi_n*pi/180; % phi_0 sorteado
        phi_0_rad = phi_0_rad + 2*pi*(F0 - F1)*tau_nB(m)/Fs; % fator de correção para F1 fora da nominal
        PHI=w1*n+w2.*(n-tau_nB(m)).*(n>=tau_nB(m))+phi_0_rad; % fase instantanea
        x=Xm.*cos(PHI);  % sinal x[n] 
        subplot(3,1,m)
        plot(n - tau_nB(m), x); hold on; %plot(tau_nB,Xm,'x')
        axis([-max(NB/2) max(NB/2) -1 1])
        xlabel('n - $\tau_{nB}$','Interpreter','latex')
        ylabel('x[n]')
        legend('N = ')
        grid on
end