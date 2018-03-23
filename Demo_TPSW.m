% Demonstration of the Two-pass Split-Window
% Last modified: 28th March, 2014
% Author: Paulo Esquef, e-mail: pesquef@lncc.br

% Uses function tpsw.m
clear all; close all; clc
% Generates an test signal: 
x=rand(1000,1); % random signal with uniform distribution
x(300:310)=x(300:310)+5; % forces a peak in x at about sample 300

% Compares the TPSW filtering with a regular moving-average filter
Nw=50; % Half the length (in samples) of the split-window
Ng=5; % Half the length (in samples) of the central gap
a=2; % 
x_tpsw=tpsw(x,Nw,Ng,a);  % the length of the split-window is 100 samples
x_ma=filter2(ones(100,1)./100,x); 

close all; plot(x); hold on; plot(x_tpsw,'r','linewidth',2); 
plot(x_ma,'k','linewidth',2);
set(gcf,'position',[277 209 1000 742]);
set(gca,'fontsize',15)
xlabel('Time (samples)')
ylabel('Magnitude')
legend('original','tpsw','moving-average')
