clear all; close all; clc;
load('boxplot_ROCOF_geral_kfest.mat')

tits = ["EF1","EF2","EF3","EF4","EF5","EF6"];
% 
% % sao 4 tipos de saltos, 2 tipos de plot (bp e hist), para 1 variaveis
% % variavel ri - 4 saltos - 2 tipos
% tipo_salto = ["salto_fase"; "salto_mag"; "salto_freq"; "senoidal"];
% var_observada = ["d_{max}"];
% 
% for i = 0:(length(var_observada)-1)
%     % variavel
%     for j=0:(length(tipo_salto)-1)
%         % tipo salto
%     filename(8*i+2*j+1) = "detector\" + var_observada(i+1) + tipo_salto(j+1)+"_plot"+ num2str(SNR)+".png";
%     filename(8*i+2*j+2) = "detector\" + var_observada(i+1) + tipo_salto(j+1)+"_boxplot"+ num2str(SNR)+".png";
%     end
% end

subs = 231:1:236;
f=0;

for i = 1:6
f = f+1;figure(f);
subplot(511)
h1 = histogram(kfest_zero(i,:), 'FaceAlpha',0.5); hold on;
h4 = histogram(kfest_freq(i,:), 'FaceAlpha',0.5); hold on;
h1.Normalization = 'probability';h1.BinWidth = 0.02;
h4.Normalization = 'probability';h4.BinWidth = 0.02;
xlim([-0.5 1.5])
legend('Senoidal','Salto Frequencia');
title("k_f estimado - EF" + i + " - SNR = "+SNR+" dB"); %ylabel('Densidade de probabilidade')

subplot(512)
h2 = histogram(kfest_fase(i,:), 'FaceAlpha',0.5); hold on;
h4 = histogram(kfest_freq(i,:), 'FaceAlpha',0.5); hold on;
h2.Normalization = 'probability';h2.BinWidth = 0.02;
h4.Normalization = 'probability';h4.BinWidth = 0.02;
xlim([-0.5 1.5])
legend('Salto Fase','Salto Frequencia');
subplot(513)
h3 = histogram(kfest_mag(i,:), 'FaceAlpha',0.5); hold on;
h4 = histogram(kfest_freq(i,:), 'FaceAlpha',0.5); hold on;
h3.Normalization = 'probability';h3.BinWidth = 0.02;
h4.Normalization = 'probability';h4.BinWidth = 0.02;
xlim([-0.5 1.5])
legend('Salto Magnitude','Salto Frequencia');
subplot(514)
h5 = histogram(kfest_fm(i,:), 'FaceAlpha',0.5); hold on;
h4 = histogram(kfest_freq(i,:), 'FaceAlpha',0.5); hold on;
h5.Normalization = 'probability';h5.BinWidth = 0.02;
h4.Normalization = 'probability';h4.BinWidth = 0.02;
xlim([-0.5 1.5])
legend('Salto FaseMagnitude','Salto Frequencia');
%title('Limiar de detecção - EF6 - SNR = 60 dB')
subplot(515)
h6 = histogram(kfest_sag(i,:), 'FaceAlpha',0.5); hold on;
h4 = histogram(kfest_freq(i,:), 'FaceAlpha',0.5); hold on;
h6.Normalization = 'probability';h6.BinWidth = 0.02;
h4.Normalization = 'probability';h4.BinWidth = 0.02;
xlim([-0.5 1.5])
legend('Afundamento','Salto Frequencia');
%title('Limiar de detecção - EF6 - SNR = 60 dB')
fname = convertStringsToChars("detector\histograma_kf_EF" + num2str(i) + "_SNR"+num2str(SNR)+".png");
saveas(gcf,fname)
end
% NAO FICOU BOM
% % separacao dos sinais em funcao de phi0
% for i = 1:6
% f = f+1;
% figure(f); plot(phi_n,kfest_sag(i,:),'b.')
% hold on; plot(phi_n,kfest_freq(i,:),'r.')
% xlabel('\phi_0 [graus]'); ylabel('k_f estimado [Hz]')
% legend('afundamento','salto frequencia')
% title("k_f estimado por EF"+i)
% end
% separacao dos sinais em funcao de tau

for i = 1:6
f = f+1;
figure(f); plot(tau_vec(1,:),kfest_sag(i,:),'b.')
hold on; plot(tau_vec(1,:),kfest_freq(i,:),'r.')
xlabel('\tau [s]'); ylabel('k_f estimado [Hz]')
legend('afundamento','salto frequencia')
title("k_f estimado por EF"+i)
end


% 
% f = f+1;figure(f); 
% boxplot(kfest_zero','Colors','kkkk'); hold on;
% boxplot(kfest_freq','Colors','bbbb');
% 
% for i=1:6
%     f = f+1;figure(f); 
%     ylabel('k_f estimado [Hz]')
%     labels = {'Sen','Fase', 'Mag', 'Freq'};
%     boxplot([kfest_zero(i,:)' kfest_fase(i,:)' kfest_mag(i,:)' kfest_freq(i,:)'],'Labels',labels);
%     title("EF" + i + " - kf estimado");
% end

fname = convertStringsToChars("detector\dados_kf_SNR"+num2str(SNR));
save(fname)