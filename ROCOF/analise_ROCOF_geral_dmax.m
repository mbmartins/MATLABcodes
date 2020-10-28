clear all; close all; clc;
load('boxplot_ROCOF_geral_dmax.mat')
%load('detector\dados_SNR80.mat')


tits = ["EF1","EF2","EF3","EF4","EF5","EF6"];

% sao 4 tipos de saltos, 2 tipos de plot (bp e hist), para 1 variaveis
% variavel ri - 4 saltos - 2 tipos
tipo_salto = ["salto_fase"; "salto_mag"; "salto_freq"; "senoidal"];
var_observada = ["d_{max}"];

for i = 0:(length(var_observada)-1)
    % variavel
    for j=0:(length(tipo_salto)-1)
        % tipo salto
    filename(8*i+2*j+1) = "detector\" + var_observada(i+1) + tipo_salto(j+1)+"_plot"+ num2str(SNR)+".png";
    filename(8*i+2*j+2) = "detector\" + var_observada(i+1) + tipo_salto(j+1)+"_boxplot"+ num2str(SNR)+".png";
    end
end

subs = 231:1:236;


f = 1; figB=figure(f); 
for k = 1:6
    subplot(subs(k)); plot(abs(riraw_fase(k,:))); %hold on; plot(ruraw_fase(k,:)); title(tits(k));
    title("EF"+k);
    xlabel('Samples'); ylabel('d[n] [Hz/s]'); %legend('ri','ru');
end
B = convertStringsToChars(filename(f));saveas(gcf,B)
figB=figure(f+1); title('Salto fase')
% for k = 1:6
%     subplot(subs(k)); histogram(riraw_fase(k,:)); title(tits(k))
% end
boxplot(draw_fase,'Labels', tits); ylabel('d_{max} [Hz/s]'); xlabel('EF'); title('Salto fase')
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2; figB=figure(f);
for k = 1:6
    subplot(subs(k)); plot(abs(riraw_mag(k,:))); %hold on; plot(ruraw_mag(k,:)); title(tits(k));
    title("EF"+k);
    xlabel('Samples'); ylabel('d[n] [Hz/s]'); %legend('ri','ru');
end
endB = convertStringsToChars(filename(f));saveas(gcf,B)%+"_boxplotFE.png")%("salto mag\salto_mag_geral"+"boxplotFE.png")
figB=figure(f+1); title('Salto magnitude')
% for k = 1:6
%     subplot(subs(k)); histogram(riraw_mag(k,:)); title(tits(k))
% end
boxplot(draw_mag,'Labels', tits); ylabel('d_{max} [Hz/s]'); xlabel('EF'); title('Salto magnitude')
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2;figB=figure(f);
for k = 1:6
    subplot(subs(k)); plot(abs(riraw_freq(k,:))); %hold on; plot(ruraw_freq(k,:)); title(tits(k));
    title("EF"+k);
    xlabel('Samples'); ylabel('d[n] [Hz/s]');%legend('ri','ru');
end
B = convertStringsToChars(filename(f));saveas(gcf,B)
figB=figure(f+1); title('Salto frequencia')
% for k = 1:6
%     subplot(subs(k)); histogram(riraw_freq(k,:)); title(tits(k))
% end
boxplot(draw_freq,'Labels', tits); ylabel('d_{max} [Hz/s]'); xlabel('EF'); title('Salto frequencia')
B = convertStringsToChars(filename(f+1));saveas(gcf,B)

f = f+2;figB=figure(f); 
for k = 1:6
    subplot(subs(k)); plot(abs(riraw_zero(k,:))); %hold on; plot(ruraw_zero(k,:));title(tits(k));
    title("EF"+k);
    xlabel('Samples'); ylabel('d[n] [Hz/s]'); %legend('ri','ru');
end
B = convertStringsToChars(filename(f));saveas(gcf,B)%+"_boxplotFE.png")%("senoidal\senoidal_geral"+"boxplotFE.png")
figB=figure(f+1); %title('Senoidal')
% for k = 1:6
%     subplot(subs(k)); histogram(riraw_zero(k,:)); title(tits(k))
% end
boxplot(draw_zero,'Labels', tits); ylabel('d_{max} [Hz/s]'); xlabel('EF'); title('Senoidal')
B = convertStringsToChars(filename(f+1));saveas(gcf,B)%+"_histFE.png")%("salto fase\salto_fase_geral"+"histFE.png")

f = f+2; figure(f); hold off;
boxplot(draw_zero,'Colors','kkkk','Symbol','kx', 'Notch','on', 'Labels', tits); hold on;
boxplot(draw_freq,'Colors','bbbb', 'Symbol','bo','Notch','on','Labels', tits);
title('Limiar para salto em frequência')
ylabel('d_{max} [Hz/s]'); 
%legend('Senoidal', 'Salto Frequencia')
fname = convertStringsToChars("detector\boxplot_comparacao_SNR"+num2str(SNR)+".png");
saveas(gcf,fname)

f = f+1;figure(f); nbins = 30;
h1 = histogram(draw_zero(:,5), 'FaceAlpha',0.5); hold on;
h2 = histogram(draw_freq(:,5), 'FaceAlpha',0.5); legend('Senoidal','Salto Frequencia')
h1.Normalization = 'probability';
h1.BinWidth = 0.02;
h2.Normalization = 'probability';
h2.BinWidth = 0.02;
title("Limiar de detecção - EF5 - SNR ="+SNR+" dB"); ylabel('Densidade de probabilidade')
fname = convertStringsToChars("detector\histograma_comparacao_EF5_SNR"+num2str(SNR)+".png");
saveas(gcf,fname)

f = f+1;figure(f); nbins = 30;
h1 = histogram(draw_zero(:,6), 'FaceAlpha',0.5); hold on;
h2 = histogram(draw_freq(:,6), 'FaceAlpha',0.5); legend('Senoidal','Salto Frequencia')
h1.Normalization = 'probability';
h1.BinWidth = 0.02;
h2.Normalization = 'probability';
h2.BinWidth = 0.02; 
title("Limiar de detecção - EF6 - SNR = "+SNR+" dB"); ylabel('Densidade de probabilidade')
%title('Limiar de detecção - EF6 - SNR = 60 dB')
fname = convertStringsToChars("detector\histograma_comparacao_EF6_SNR"+num2str(SNR)+".png");
saveas(gcf,fname)

f = f+1;figure(f); nbins = 30;
h1 = histogram(draw_zero(:,6), 'FaceAlpha',0.5); hold on;
h2 = histogram(draw_freq(:,6), 'FaceAlpha',0.5); legend('Senoidal','Salto Frequencia')
h1.Normalization = 'probability';
h1.BinWidth = 0.02;
h2.Normalization = 'probability';
h2.BinWidth = 0.02; 
title("Limiar de detecção - EF6 - SNR = "+SNR+" dB"); ylabel('Densidade de probabilidade')
%title('Limiar de detecção - EF6 - SNR = 60 dB')
fname = convertStringsToChars("detector\histograma_comparacao_EF6_SNR"+num2str(SNR)+".png");
saveas(gcf,fname)

% boxplot(kfest_zero','Colors','kkkk'); hold on;
% boxplot(kfest_freq','Colors','bbbb'); hold on;


fname = convertStringsToChars("detector\dados_SNR"+num2str(SNR));
save(fname)