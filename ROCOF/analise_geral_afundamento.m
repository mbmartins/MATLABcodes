clear all; close all; clc;
load('afundamento_geral.mat')

tits = ["EF1","EF2","EF3","EF4","EF5","EF6"];

% sao 1 tipos de saltos, 2 tipos de plot (bp e hist), para 3 variaveis
variavel = ["FE"; "FE1"; "KFE"];

for i = 0:2 
    % variavel
    filename(2*i+1) = "geral\afundamento" + variavel(i+1) + "_boxplot.png";
    filename(2*i+2) = "geral\afundamento" + variavel(i+1) + "_hist.png";
end

subs = 231:1:236;

f = 1; figB=figure(f); boxplot(FEraw_fase')
title('Afundamento');xlabel('EF'); ylabel('FE [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%+"_boxplotFE.png")%("Afundamento\salto_fase_geral"+"boxplotFE.png")
figB=figure(f+1); title('Afundamento')
for k = 1:6
    subplot(subs(k)); histogram(FEraw_fase(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("Afundamento\salto_fase_geral"+"histFE.png")

% ---- figuras para FE1
f = f+2;figB=figure(f); boxplot(fE1raw_fase')
title('Afundamento');xlabel('EF'); ylabel('FE_1 [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%)%(filename(f))
figB=figure(f+1); title('Afundamento')
for k = 1:6
    subplot(subs(k)); histogram(fE1raw_fase(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("Afundamento\salto_fase_geral"+"histFE.png")

%--- figuras para KFE-----------------------
f = f+2;figB=figure(f); boxplot(fE1raw_fase')
title('Afundamento');xlabel('EF'); ylabel('K_fE [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%)%(filename(f))
figB=figure(f+1); title('Afundamento')
for k = 1:6
    subplot(subs(k)); histogram(kfEraw_fase(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("Afundamento\salto_fase_geral"+"histFE.png")

