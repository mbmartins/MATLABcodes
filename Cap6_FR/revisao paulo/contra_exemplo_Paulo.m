%Por que n�o? Faz sentido as distribui��es de k_f para saltos de magnitude que voc� encontrou?

%Minha expectativa era que a distribui��o de k_f para salto de magnitude sofresse um bias na m�dia, mas fosse concentrada. 

%Verificar os experimentos.  Para mim, no caso de salto de magnitude, n�o faz sentido usar a compensa��o pela magnitude a_i[n] no c�lculo de f_i[n], como em (5.57) e (5.58). Cria-se um distor��o em degrau artificial na fun��o de frequ�ncia instant�nea que torna os resultados sem sentido.  

%Explica a diferen�a entre os resultados obtidos pela simula��o abaixo e os que voc� encontrou! 

%Contra-exemplo que mostra a separa��o entre as distribui��es de k_f s�o separadas, para estimadores de f_1 e f_2, similares aos de (5.57) e (5.58) e \hat{k}_f=\hat{f}_2 - \hat{f}_2, s� que no caso da estimativa de k_f para o caso de salto de magnitude, a compensa��o por a_i[n] n�o � obviamente usada. 


clear all; close all; clc
fs=5000;% frequencia de amostragem 
f1=60; % frequencia nominal em Hz
kf=1;  % salto em frequ�ncia em Hz
kx=0;  % em percentual de A (salto em mag)
ka=0;  % em graus (salto em fase em graus)
A=1; % magnitude nominal 
fasei=0; % fase inicial em graus
phi=fasei*pi/180; % fase inicial em rad 
N=500; 
tau=250;  % tau_n
n=(0:(N-1))'; % base de tempo unit�rio
w1=f1/fs;
wkf=kf/fs;
fase=2*pi*w1.*n+phi; % fase limpa
mag=A.*ones(N,1);  % magnitude limpa

if abs(kf),
    fase_extra=2*pi*wkf.*(n-tau);
    fase(tau+1:end)=fase(tau+1:end)+fase_extra(tau+1:end);
end

if abs(ka),
    fase(tau:end)=fase(tau:end)+ka*pi/180;
end

if abs(kx),
    mag(tau:end)=mag(tau:end)+A*kx;
end
   
SNRdB=60; % SRN alvo em dB
x=mag.*cos(fase); % sinal sem ru�do

vx=var(x); % vari�ncia do ru�do
q=round(0.05*N);
MC=1000;
eka=zeros(MC,1);
ekf=zeros(MC,1);
ekx=zeros(MC,1);
g=fs/(2*pi);
for jj=1:MC,
cn=sqrt(vx/10^(SNRdB/10))*randn(size(x)); % componente de ru�do
xn=x+cn; % sinal com ru�do

xa=hilbert(xn); 
%psi=phase(xa);
psi=unwrap(angle(xa));
ai=abs(xa);
fi=gradient(psi);  % fi do sinal ruidoso 
if ~kf,
    aic=fi.*ai./median(fi);
elseif kf, 
    aic=ai;
end

if abs(kx)==0,
    fic=ai.*fi./median(ai);
elseif abs(kx)>=0,
    fic=fi;
end


ef1=median(fic(q:tau-q));
ef2=median(fic(tau+q:end));
ekf(jj)=(ef2-ef1)*g;

if 0
epsic=psi-ef1.*n-phi;
epsi1=median(epsic(q:tau-q));
epsi2=median(epsic(tau+q:end));
eka(jj)=(epsi2-epsi1)*180/pi;
end

%plot(fic.*g)
%pause
if 1
mfic=median(fic(q:end-q));
aux=max(abs(fic(q:end-q)-mfic));
eka(jj)=aux*g/(2*pi);
end

ea1=median(aic(q:tau-q));
ea2=median(aic(tau+q:end));
ekx(jj)=ea2-ea1;
end

ekf_freq=ekf;

%%%%%%%%%%%%%%%%%%%%%

kf=0;  % salto em frequ�ncia em Hz
kx=1.1;  % em percentual de A (salto em mag)
ka=0;  % em graus (salto em fase em graus)
A=1; % magnitude nominal 
fasei=0; % fase inicial em graus
phi=fasei*pi/180; % fase inicial em rad 
N=500; 
%pertau=.5;  % localiza��o do salto em termos do percentual do tamanho da janela 
%tau=round(N*pertau); % �ndice do instante (arredondado) de ocorr�ncia do salto
tau=250;  % tau ajustado para o salto de fase coincidir com um maximo local do seno
n=(0:(N-1))'; % base de tempo unit�rio
w1=f1/fs;
wkf=kf/fs;
fase=2*pi*w1.*n+phi; % fase limpa
mag=A.*ones(N,1);  % magnitude limpa

if abs(kf),
    fase_extra=2*pi*wkf.*(n-tau);
    fase(tau+1:end)=fase(tau+1:end)+fase_extra(tau+1:end);
end

if abs(ka),
    fase(tau:end)=fase(tau:end)+ka*pi/180;
end

if abs(kx),
    mag(tau:end)=mag(tau:end)+A*kx;
end
   
SNRdB=60; % SRN alvo em dB
x=mag.*cos(fase); % sinal sem ru�do

vx=var(x); % vari�ncia do ru�do
q=round(0.05*N);
MC=1000;
eka=zeros(MC,1);
ekf=zeros(MC,1);
ekx=zeros(MC,1);
g=fs/(2*pi);
for jj=1:MC,
cn=sqrt(vx/10^(SNRdB/10))*randn(size(x)); % componente de ru�do
xn=x+cn; % sinal com ru�do

xa=hilbert(xn); 
%psi=phase(xa);
psi=unwrap(angle(xa));
ai=abs(xa);
fi=gradient(psi);  % fi do sinal ruidoso 
if ~kf,
    aic=fi.*ai./median(fi);
elseif kf, 
    aic=ai;
end

if abs(kx)==0,
    fic=ai.*fi./median(ai);
elseif abs(kx)>=0,
    fic=fi;
end

ef1=median(fic(q:tau-q));
ef2=median(fic(tau+q:end));
ekf(jj)=(ef2-ef1)*g;


if 0
epsic=psi-ef1.*n-phi;
epsi1=median(epsic(q:tau-q));
epsi2=median(epsic(tau+q:end));
eka(jj)=(epsi2-epsi1)*180/pi;
end


if 1
mfic=median(fic(q:end-q));
aux=max(abs(fic(q:end-q)-mfic));
eka(jj)=aux*g/(2*pi);
end

ea1=median(aic(q:tau-q));
ea2=median(aic(tau+q:end));
ekx(jj)=ea2-ea1;

end

ekf_mag=ekf;
hist(ekf_freq); hold on;  hist(ekf_mag);




