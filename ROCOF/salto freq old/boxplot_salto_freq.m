clear all; close all; clc;
load('boxplot_ROCOF_geral_dmax.mat')
%load('detector\dados_SNR40.mat')

tits = ["EF1";"EF2";"EF3";"EF4";"EF5";"EF6"];
subs = 321:1:326;

figs(1) = figure('Units','normalized','Position',[0 0 0.5 1]);
for i = 1:6
    draw_ALL = [draw_zero(:,i) draw_freq(:,i)];
    
    subplot(subs(i));
    boxplot(draw_ALL, 'Labels', {'Senoidal', 'Salto Frequência'});
    ylabel('d_{rmax} [Hz/s]'); xlabel(tits(i)) 
    grid on;
end
    fname = convertStringsToChars("boxplot_seno_"+num2str(SNR)+".png");
    saveas(gcf,fname)

f = i+1;figure(f); nbins = 30;
h1 = histogram(draw_zero(:,4), 'FaceAlpha',0.5); hold on;
h2 = histogram(draw_freq(:,4), 'FaceAlpha',0.5); legend('Senoidal','Salto Frequencia')
h1.Normalization = 'probability';
h1.BinWidth = 1.e-4;
h2.Normalization = 'probability';
h2.BinWidth = 1.e-4;
title("Limiar de detecção - EF4 - SNR ="+SNR+" dB"); ylabel('Densidade de probabilidade')
fname = convertStringsToChars("histograma_comparacao_EF4_SNR"+num2str(SNR)+".png");
saveas(gcf,fname)    
    
      
figs(i+6) = figure('Units','normalized','Position',[0 0 0.7 0.8]);    
for i = 1:6
    draw_ALL = [draw_freq(:,i) draw_mag(:,i) draw_fase(:,i) draw_sag(:,i)];
    subplot(subs(i));
    boxplot(draw_ALL, 'Labels', {'Frequência', 'Magnitude', 'Fase', 'Afundamento'});
    ylabel('d_{rmax} [Hz/s]'); xlabel(tits(i)) 
    grid on;
end

fname = convertStringsToChars("boxplotALLoutros"+num2str(SNR)+".png");
saveas(gcf,fname)
