%testes PATV
close all; clear all; clc;

t = 0:10;
u_d = [0 0 0 0 0 1 1 1 1 1 1];
a = 1;
b = 1;
y = a*t + b + u_d;
Nit = 30;
lambda = 1e-6; 
%lambda pequeno significa que o degrau (sinal x) � 'extra�do' de forma mais
%pronunciada do sinal de entrada 'y'. Significa que os par�metros do
%polin�mio p estimado s�o mais pr�ximos dos par�metros da parcela
%polinomial do sinal y. Provavelmente, lambda menor d� uma estimativa
%melhor da inclina��o, o que pode ser �til para estimar a frequ�ncia.

%todavia, nos sinais provenientes do sinal anal�tico para o nosso
%experimento, n�o h� somente ru�do, mas oscila��es provenientes das n�o
%idealidades. quando reduzimos muito o lambda, o sinal x tenta incorporar estas
%oscila��es e o ru�do. O valor m�dio antes de tau e depois de n�o � semelhante, 
%o que afeta a inclina��o da estimativa p de theta.

d = 1;

[x, p, cost, u, v] = patv_MM(y, d, lambda, Nit);
plot(cost); title('Cost Function')
figure
plot(t,p,t,x,t,y)
legend('p','x','y')

a_est = p(2)-p(1)
erro = a - a_est

% conclus�o: o fator constante � absorvido em p, n�o em x. Isso � uma
% diferen�a da implementa��o do artigo de 2013, em rela��o ao de 2012.
