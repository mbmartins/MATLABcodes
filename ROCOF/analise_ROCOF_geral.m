clear all; close all; clc;
load('boxplot_ROCOF_geral.mat')

tits = ["EF1","EF2","EF3","EF4","EF5","EF6"];

% sao 4 tipos de saltos, 2 tipos de plot (bp e hist), para 1 variaveis
% variavel ri - 4 saltos - 2 tipos
tipo_salto = ["salto_fase"; "salto_mag"; "salto_freq"; "senoidal"];
var_observada = ["ri", "d_{max}"];

for i = 0:(length(var_observada)-1)
    % variavel
    for j=0:(length(tipo_salto)-1)
        % tipo salto
    filename(8*i+2*j+1) = "geral\" + var_observada(i+1) + tipo_salto(j+1)+"_plot.png";
    filename(8*i+2*j+2) = "geral\" + var_observada(i+1) + tipo_salto(j+1)+"_boxplot.png";
    end
end

subs = 231:1:236;

f = 1; figB=figure(f); 
for k = 1:6
    subplot(subs(k)); plot(riraw_fase(k,:)); %hold on; plot(ruraw_fase(k,:)); title(tits(k));
    title("EF"+k);
    xlabel('Samples'); ylabel('r_i [Hz/s]'); %legend('ri','ru');
end
B = convertStringsToChars(filename(f));saveas(gcf,B)
figB=figure(f+1); title('Salto fase')
% for k = 1:6
%     subplot(subs(k)); histogram(riraw_fase(k,:)); title(tits(k))
% end
boxplot(draw_fase); ylabel('r_i [Hz/s]'); xlabel('EF')
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2; figB=figure(f);
for k = 1:6
    subplot(subs(k)); plot(riraw_mag(k,:)); %hold on; plot(ruraw_mag(k,:)); title(tits(k));
        title("EF"+k);
    xlabel('Samples'); ylabel('r_i [Hz/s]'); %legend('ri','ru');
end
endB = convertStringsToChars(filename(f));saveas(gcf,B)%+"_boxplotFE.png")%("salto mag\salto_mag_geral"+"boxplotFE.png")
figB=figure(f+1); title('Salto magnitude')
% for k = 1:6
%     subplot(subs(k)); histogram(riraw_mag(k,:)); title(tits(k))
% end
boxplot(draw_mag); ylabel('r_i [Hz/s]'); xlabel('EF')
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2;figB=figure(f);
for k = 1:6
    subplot(subs(k)); plot(riraw_freq(k,:)); %hold on; plot(ruraw_freq(k,:)); title(tits(k));
        title("EF"+k);
    xlabel('Samples'); ylabel('r_i [Hz/s]');%legend('ri','ru');
end
B = convertStringsToChars(filename(f));saveas(gcf,B)
figB=figure(f+1); title('Salto frequencia')
% for k = 1:6
%     subplot(subs(k)); histogram(riraw_freq(k,:)); title(tits(k))
% end
boxplot(draw_freq); ylabel('r_i [Hz/s]'); xlabel('EF')
B = convertStringsToChars(filename(f+1));saveas(gcf,B)

f = f+2;figB=figure(f); 
for k = 1:6
    subplot(subs(k)); plot(riraw_zero(k,:)); %hold on; plot(ruraw_zero(k,:));title(tits(k));
        title("EF"+k);
    xlabel('Samples'); ylabel('r_i [Hz/s]'); %legend('ri','ru');
end
B = convertStringsToChars(filename(f));saveas(gcf,B)%+"_boxplotFE.png")%("senoidal\senoidal_geral"+"boxplotFE.png")
figB=figure(f+1); %title('Senoidal')
% for k = 1:6
%     subplot(subs(k)); histogram(riraw_zero(k,:)); title(tits(k))
% end
boxplot(draw_zero); ylabel('r_i [Hz/s]'); xlabel('EF')
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

