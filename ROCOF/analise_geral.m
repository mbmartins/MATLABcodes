clear all; close all; clc;
load('boxplot_geral.mat')

tits = ["EF1","EF2","EF3","EF4","EF5","EF6"];

% sao 4 tipos de saltos, 2 tipos de plot (bp e hist), para 3 variaveis
% variavel FE - 4 saltos - 2 tipos
tipo_salto = ["salto_fase"; "salto_mag"; "salto_freq"; "senoidal"];
variavel = ["FE"; "FE1"; "KFE"];

for i = 0:2 
    % variavel
    for j=0:3 
        % tipo salto
    filename(8*i+2*j+1) = "geral\" + variavel(i+1) + tipo_salto(j+1)+"_boxplot.png";
    filename(8*i+2*j+2) = "geral\" + variavel(i+1) + tipo_salto(j+1)+"_hist.png";
    end
end


subs = 231:1:236;

f = 1; figB=figure(f); boxplot(FEraw_fase')
title('Salto fase');xlabel('EF'); ylabel('FE [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%+"_boxplotFE.png")%("salto fase\salto_fase_geral"+"boxplotFE.png")
figB=figure(f+1); title('Salto fase')
for k = 1:6
    subplot(subs(k)); histogram(FEraw_fase(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2; figB=figure(f);boxplot(FEraw_mag')
title('Salto magnitude');xlabel('EF');ylabel('FE [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%+"_boxplotFE.png")%("salto mag\salto_mag_geral"+"boxplotFE.png")
figB=figure(f+1); title('Salto magnitude')
for k = 1:6
    subplot(subs(k)); histogram(FEraw_mag(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2;figB=figure(f);boxplot(FEraw_freq')
title('Salto frequencia');xlabel('EF');ylabel('FE [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%+"_boxplotFE.png")%("salto freq\salto_freq_geral"+"boxplotFE.png")
figB=figure(f+1); title('Salto frequencia')
for k = 1:6
    subplot(subs(k)); histogram(FEraw_freq(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2;figB=figure(f); boxplot(FEraw_zero')
title('Senoidal');xlabel('EF'); ylabel('FE [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%+"_boxplotFE.png")%("senoidal\senoidal_geral"+"boxplotFE.png")
figB=figure(f+1); %title('Senoidal')
for k = 1:6
    subplot(subs(k)); histogram(FEraw_zero(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

% ---- figuras para FE1
f = f+2;figB=figure(f); boxplot(fE1raw_fase')
title('Salto fase');xlabel('EF'); ylabel('FE_1 [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%)%(filename(f))
figB=figure(f+1); title('Salto fase')
for k = 1:6
    subplot(subs(k)); histogram(fE1raw_fase(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2;figB=figure(f);boxplot(fE1raw_mag')
title('Salto magnitude');xlabel('EF');ylabel('FE_1 [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%)%("salto mag\salto_mag_geral"+"boxplotFE1.png")
figB=figure(f+1); title('Salto magnitude')
for k = 1:6
    subplot(subs(k)); histogram(fE1raw_mag(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2; figB=figure(f);boxplot(fE1raw_freq')
title('Salto frequencia');xlabel('EF');ylabel('FE_1 [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%)%("salto_freq\salto_freq_geral"+"boxplotFE1.png")
figB=figure(f+1); title('Salto frequencia')
for k = 1:6
    subplot(subs(k)); histogram(fE1raw_freq(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2;figB=figure(f); boxplot(fE1raw_zero')
title('Senoidal');xlabel('EF'); ylabel('FE_1 [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%)%("senoidal\senoidal_geral"+"boxplotFE1.png")
figB=figure(f+1); title('Salto fase')
for k = 1:6
    subplot(subs(k)); histogram(FEraw_zero(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

%--- figuras para KFE-----------------------
f = f+2;figB=figure(f); boxplot(fE1raw_fase')
title('Salto fase');xlabel('EF'); ylabel('K_fE [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%)%(filename(f))
figB=figure(f+1); title('Salto fase')
for k = 1:6
    subplot(subs(k)); histogram(kfEraw_fase(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2;figB=figure(f);boxplot(kfEraw_mag')
title('Salto magnitude');xlabel('EF');ylabel('K_fE [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%)%("salto mag\salto_mag_geral"+"boxplotkfE.png")
figB=figure(f+1); title('Salto magnitude')
for k = 1:6
    subplot(subs(k)); histogram(kfEraw_mag(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2; figB=figure(f);boxplot(kfEraw_freq')
title('Salto frequencia');xlabel('EF');ylabel('K_fE [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%)%("salto_freq\salto_freq_geral"+"boxplotkfE.png")
figB=figure(f+1); title('Salto frequencia')
for k = 1:6
    subplot(subs(k)); histogram(kfEraw_freq(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2;figB=figure(f); boxplot(kfEraw_zero')
title('Senoidal');xlabel('EF'); ylabel('K_fE [Hz]')
B = convertStringsToChars(filename(f));saveas(gcf,B)%)%("senoidal\senoidal_geral"+"boxplotkfE.png")
figB=figure(f+1); title('Salto fase')
for k = 1:6
    subplot(subs(k)); histogram(kfEraw_zero(k,:)); title(tits(k))
end
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")
