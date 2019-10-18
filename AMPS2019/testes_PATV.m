%testes PATV
close all; clear all; clc;

t = 0:10;
u_d = [0 0 0 0 0 1 1 1 1 1 1];
a = 1;
b = 1;
y = a*t + b + u_d;
Nit = 30;
lambda = 1e-6; 
%lambda pequeno significa que o degrau (sinal x) é 'extraído' de forma mais
%pronunciada do sinal de entrada 'y'. Significa que os parâmetros do
%polinômio p estimado são mais próximos dos parâmetros da parcela
%polinomial do sinal y. Provavelmente, lambda menor dá uma estimativa
%melhor da inclinação, o que pode ser útil para estimar a frequência.

%todavia, nos sinais provenientes do sinal analítico para o nosso
%experimento, não há somente ruído, mas oscilações provenientes das não
%idealidades. quando reduzimos muito o lambda, o sinal x tenta incorporar estas
%oscilações e o ruído. O valor médio antes de tau e depois de não é semelhante, 
%o que afeta a inclinação da estimativa p de theta.

d = 1;

[x, p, cost, u, v] = patv_MM(y, d, lambda, Nit);
plot(cost); title('Cost Function')
figure
plot(t,p,t,x,t,y)
legend('p','x','y')

a_est = p(2)-p(1)
erro = a - a_est

% conclusão: o fator constante é absorvido em p, não em x. Isso é uma
% diferença da implementação do artigo de 2013, em relação ao de 2012.
