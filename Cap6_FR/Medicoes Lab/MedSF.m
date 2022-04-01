function [f1,f2,f_r,h_f] = MedSF(Signal, tau_n,Fs)
%EF1
% 1 - amostragem ideal
% 2 - Psi obtido pela fase do sinal analitico
% 3 - f_u[n] obtida pelo gradient de Psi
% 4 - considera f_u = f_i
% 5 
%   - retira 5% das amostras proximas das bordas
%   - considera f_r = media de f_u, aproximada pela mediana
%   - estima f_1 = mediana de f_u até tau_n
%   - estima f_2 = mediana de f_u a partir de tau_n
% 6 - estima r_i[n] pelo gradiente de f_i[n]

    z=hilbert(Signal');
    Psi_i = unwrap(angle(z)); %[rad]
    f_i=gradient(Psi_i)*Fs/(2*pi);% Hilbert estimate of the instantaneous frequency of z
    NSamples = length(Signal);
    br = 0.05; %80/480; % fraction of samples to be ignored 
    brn = floor(br*NSamples); % number of samples to be ignored
    brmask1 = [(brn+1:tau_n)];
    brmask2 = [(tau_n+1:NSamples-brn)];
    %brmask = [brmask1 brmask2];
    %f_u_trunc = f_i(brmask);
    %figure; plot(f_u_trunc); title('Frequencia instantanea')
    
    %ri = zeros(NSamples,1);
    %ri(brmask) = gradient(f_u_trunc);
    f1 = median(f_i(brmask1));
    f2 = median(f_i(brmask2));
    tau = tau_n/NSamples;
    f_r = tau*f1 + (1-tau)*f2;
    h_f = f2 - f1;
    
% ---- Debug ----
% subplot(1,2,1)
% plot(f_u_trunc)
% xlabel('Samples'); ylabel('f_u [Hz]')
% subplot(1,2,2)
% h = histogram(f_u_trunc,'Orientation','horizontal');
% h.NumBins = 30;
% xlabel('Occurrences'); ylabel('f_u [Hz]')

%     figure;
%     nn = brn:NSamples-brn-1;
%     plot(nn,ri);
%     xlabel('x','Interpreter','latex');
%     ylabel('y','Interpreter','latex');
%     xlabel('Samples')
%     ylabel('$\hat{r}_i[n]$ [Hz/s]')

%coeficiente de assimetria
%sk = skewness(f_u_trunc,0);

