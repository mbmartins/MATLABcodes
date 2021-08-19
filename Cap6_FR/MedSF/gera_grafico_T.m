clear all; close all; clc;

%geracao dos graficos
c = ['k','b','c','r','m','g'];

plot_type = 6;
fsize = 16; %font sizes

pasta = "figuras_hf_1_10000\";
load(pasta+"estimadoresFR_salto_freq_T.mat");
VEC = Ncycles_vec/60;
sxlabel = 'T [s]';
var = "T";
vxlim = [3 30]/60;
KE = KFE;

figure(plot_type)
figALL(plot_type) = figure(plot_type);
figALL(plot_type).Units = 'normalized';
figALL(plot_type).Position = [0 0 1.0 1.0];
% ---- plotar graficos
k = 1; %MedSF = antigo EF1
    for i =1:2
        subplot(2,2,1)
        plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
        set(gca,'FontSize',fsize);
        pleg(k) = plot(VEC,FE(:,k),c(k)+"-"); hold on;
        shade(VEC,FE(:,k) - FE_std(:,k),[c(k),'--'],VEC,FE(:,k) + FE_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
        xlabel(sxlabel,'FontSize',fsize); ylabel('FE [Hz]','FontSize',fsize)
        xlim(vxlim); %ylim([-0.06 0.06])
        grid on; grid minor;
        subplot(2,2,2)
        plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
        set(gca,'FontSize',fsize);
        pleg(k+1) = plot(VEC,KE(:,k),c(k)+"-"); hold on;
        shade(VEC,KE(:,k) - kf_std(:,k),[c(k),'-.'],VEC,KE(:,k) + kf_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
        xlabel(sxlabel,'FontSize',fsize); ylabel('h_fE [Hz]','FontSize',fsize)
        xlim(vxlim); %ylim([-1.1 0.2])
        grid on; grid minor;
        subplot(2,2,3)
        plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
        set(gca,'FontSize',fsize);
        pleg(k+2) = plot(VEC,FE1(:,k),c(k)+"-"); hold on;
        shade(VEC,FE1(:,k) - f1_std(:,k),[c(k),'-.'],VEC,FE1(:,k) + f1_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
        xlabel(sxlabel,'FontSize',fsize); ylabel('FE_1 [Hz]','FontSize',fsize)
        xlim(vxlim); %ylim([-0.1 0.5])
        grid on; grid minor;
        %k=k+4; %MedSF-PATV
        subplot(2,2,4)
        plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
        set(gca,'FontSize',fsize);
        pleg(k+3) = plot(VEC,FE2(:,k),c(k)+"-"); hold on;
        shade(VEC,FE2(:,k) - f2_std(:,k),[c(k),'-.'],VEC,FE2(:,k) + f2_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
        xlabel(sxlabel,'FontSize',fsize); ylabel('FE_2 [Hz]','FontSize',fsize)
        xlim(vxlim); %ylim([-0.1 0.5])
        grid on; grid minor;
        
        k=k+4; %next loop: MedSF-PATV antigo EF5, por isso k+4
        
    end
    k = 1;
    for i = 1:4
        lgd(i) = legend([pleg(k) pleg(k+4)],"MedSF","MedSF-PATV");
        lgd(i).Location='south';
        lgd(i).Orientation='horizontal';
        lgd(i).Box='off';
        %legend('boxoff')
        lgd(i).FontSize = fsize;
        k = k+1;
    end

% --- salvar figuras
savefig(figALL(plot_type),pasta+"salto_freq_ALL_"+var+"_MedSF")
saveas(figALL(plot_type),pasta+"salto_freq_ALL_"+var+"_MedSF.png")

message = plot_type + " graficos criados na pasta "+ pasta;