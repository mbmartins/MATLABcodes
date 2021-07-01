%geracao dos graficos
clc; clear all; close all; c = ['k','b','c','r','m','g'];

% casos:
% 1 - phi_0
% 2 - tau
% 3 - f1
% 4 - altura do salto hf
% 5 - frequencia de amostragem fs
% 6 - tamanho da janela T

plot_type = 6; 
fsize = 16;
    
switch plot_type
    case 1
    % ---- Fase inicial -----  influencia de phi0
    load('salto_freq_phi0.mat'); close all;
    var="phi0"; VEC = phi_n;
    sxlabel = '\phi_0 [graus]';
    vxlim = [0 175];
    FE = FE_ruido;
    KE = KE_ruido;
    FE1 = FE_ruido;
    case 2
        % --- tau
    load('estimadoresFR_salto_freq_tau_n.mat');close all;  
    var="tau"; VEC=tau_n;
    sxlabel = '\tau_n [amostras]';
    vxlim = [48 432];
    KE = KFE1;
    case 3
        % ---- f1
    load('estimadoresFR_salto_freq_F1.mat');close all;
    VEC = F1_vec;
    var = "F1"
    sxlabel = 'f_1 [Hz]';
    vxlim = [58 62];
    KE = KFE;
    case 4
        % ---- altura do salto hf
        load('estimadoresFR_salto_freq_kf.mat');close all;
        VEC = kf_vec;
        sxlabel = 'h_f [Hz]';
        var = "hf"
        vxlim = [0 3];
        KE = KFE;
        
    case 5
        % ---- frequencia de amostragem fs
        load('estimadoresFR_salto_freq_Fs.mat');close all;
        VEC = Fs_vec;
        sxlabel = 'f_s [Hz]';
        var = "fs";
        vxlim = [1000 20000];
        KE = KFE;
    case 6 
        % ---- tamanho da janela T
        load('estimadoresFR_salto_freq_T.mat');close all;
        VEC = Ncycles_vec;
        sxlabel = 'T [ciclos]';
        var = "T";
        vxlim = [3 30];
        KE = KFE;
end

figALL = figure('Units','normalized','Position',[0 0 1.0 1.0]);
% ---- plotar graficos
k = 1; %MedSF = antigo EF1
    for i =1:2
        subplot(1,3,1)
        plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
        set(gca,'FontSize',fsize);
        pleg(k) = plot(VEC,FE(:,k),c(k)+"-"); hold on;
        shade(VEC,FE(:,k) - FE_std(:,k),[c(k),'--'],VEC,FE(:,k) + FE_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
        xlabel(sxlabel,'FontSize',fsize); ylabel('FE [Hz]','FontSize',fsize)
        xlim(vxlim); %ylim([-0.06 0.06])
        grid on; grid minor;
        subplot(1,3,2)
        plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
        set(gca,'FontSize',fsize);
        pleg(k+1) = plot(VEC,KE(:,k),c(k)+"-"); hold on;
        shade(VEC,KE(:,k) - kf_std(:,k),[c(k),'-.'],VEC,KE(:,k) + kf_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
        xlabel(sxlabel,'FontSize',fsize); ylabel('h_fE [Hz]','FontSize',fsize)
        xlim(vxlim); %ylim([-1.1 0.2])
        grid on; grid minor;
        subplot(1,3,3)
        plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
        set(gca,'FontSize',fsize);
        pleg(k+2) = plot(VEC,FE1(:,k),c(k)+"-"); hold on;
        shade(VEC,FE1(:,k) - f1_std(:,k),[c(k),'-.'],VEC,FE1(:,k) + f1_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
        xlabel(sxlabel,'FontSize',fsize); ylabel('FE_1 [Hz]','FontSize',fsize)
        xlim(vxlim); %ylim([-0.1 0.5])
        grid on; grid minor;
        k=k+4; %MedSF-PATV
        
    end
    k = 1;
    for i = 1:3
        lgd(i) = legend([pleg(k) pleg(k+4)],"MedSF","MedSF-PATV");
        lgd(i).Location='southoutside';
        lgd(i).Orientation='horizontal';
        %legend('boxoff')
        lgd(i).FontSize = fsize;
        k = k+1;
    end

% --- salvar figuras
savefig(figALL,"MedSF\salto_freq_ALL_"+var+"_MedSF")
saveas(figALL,"MedSF\salto_freq_ALL_"+var+"_MedSF.png")

