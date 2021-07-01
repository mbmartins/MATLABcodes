% % gera dados
% clear all; clc; close all;
% MCiter = 10000;
% h_f =1;
% pasta = gera_dados_MC(MCiter,h_f)
% %gera graficos
% message = graficos_ALL_MedSF(pasta)

% gera dados
clear all; clc; close all;
MCiter = 10000;
h_f =-1;
pasta = gera_dados_MC(MCiter,h_f)
%gera graficos
message = graficos_ALL_MedSF(pasta)

%dados robustez
clear all; close all; clc;
MCiter = 10000;
h_f = 1;
phi_n = 180*rand(1,MCiter); %distribuicao de phi_0
tau_vec = 0.1 + (0.8)*rand(1,MCiter); %distribuicao de tau
filename = Robustez_MedSF(phi_n,tau_vec,MCiter,h_f);
message = analise_robustez(filename)

% para salto negativo
clear all; close all; clc;
MCiter = 10000;
phi_n = 180*rand(1,MCiter); %distribuicao de phi_0
tau_vec = 0.1 + (0.8)*rand(1,MCiter); %distribuicao de tau
h_f = -1;
filename = Robustez_MedSF(phi_n,tau_vec,MCiter,h_f);
message = analise_robustez(filename)

% se quiser rodar para valores fixos
%tau fixo em 0.5T
%phi_0 fixo em zero:
% phi_n = 0*ones(1,MCiter);
% tau_vec = 0.5*ones(1,MCiter);
