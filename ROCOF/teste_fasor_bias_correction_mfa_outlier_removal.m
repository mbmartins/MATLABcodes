% teste de compensação de fi pela magnitude instantânea
% uso de recursão de remoção de outlier

clear all; close all; clc

N=600; % tamanho da janela em amostras
n=0:N-1; n=n(:); % base de tempo
w0=3*pi/125; % freq fundamental
h=10*pi/180;  % salto na fase em graus
s=h.*(n>=210+20); 
fase=w0*n+s;
fase_ref=(1.0*w0)*n;
fasor=cos(fase);
vfasor=var(fasor);
SNR=60;  
vruido=vfasor./(10^(SNR/10));
fasor=fasor+sqrt(vruido)*randn(size(n));

fa=hilbert(fasor); % fasor analítico
%fi=(phase(fa));  % fase instantânea
fi=unwrap(angle(fa));  % fase instantânea
freqi=gradient(fi);
freqi=freqi(20:end-20); s=s(20:end-20); % cropped
freqi_orig=freqi;
mfa=abs(fa); 
mfa=mfa(20:end-20);
freqi=freqi.*(mfa./mean(mfa));

plot(freqi,'r')
hold on;
plot(freqi_orig,'k')
xlabel('amostras')
ylabel('amplitude')
% recursão não-linear de remoção de outlier
freqi_proc=zeros(size(freqi));
freqi_proc(1)=median(freqi);
lambda=1.0075;
for k=2:length(freqi),
    freqi_proc(k)=freqi_proc(k-1).*lambda^sign(freqi(k)-freqi_proc(k-1));
end

plot(freqi_proc,'m','linewidth',3)
legend('fi compensada magnitude','fi original','fi outlier removal')


([median(freqi_orig) median(freqi) median(freqi_proc)]-w0)*1000./w0

% para 60 dB a estimativa de freq sem qualquer compensação e/ou
% processamento apresenta melhor resultado, em termos de menor polarização

% outra maneira de filtrar o degrau na fase
phase_i = cumtrapz(freqi_orig);
figure; plot(phase_i)
phase_i_comp = cumtrapz(freqi_proc);
hold on; plot(phase_i_comp)
