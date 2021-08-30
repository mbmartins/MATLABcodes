clear all; close all; clc;

eps = [-1*ones(1,210), 0*ones(1,450), 1*ones(1,250), 2*ones(1,40)];

h = histogram(eps);
h.BinLimits = [-2 2];
h.Normalization = 'probability'
h.NumBins = 10;
ylabel('Frequência relativa','FontSize',14)
xlabel('$\epsilon_\tau$ [$\Delta t$]','Interpreter','latex','FontSize',16)
