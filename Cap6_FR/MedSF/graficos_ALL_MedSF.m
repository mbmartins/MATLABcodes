function message = graficos_ALL_MedSF(pasta)
%geracao dos graficos
c = ['k','b','c','r','m','g'];

%plot_type = 6; 
fsize = 16; %font sizes

%for plot_type = 1:7
for plot_type = 3
% casos:
% 1 - phi_0
% 2 - tau
% 3 - f1
% 4 - altura do salto hf
% 5 - frequencia de amostragem fs
% 6 - tamanho da janela T
% 7 - erro de estimativa de tau

switch plot_type
    case 1
    % ---- Fase inicial -----  influencia de phi0
    load(pasta+"salto_freq_phi0.mat"); 
    var="phi0"; VEC = phi_n;
    sxlabel = '\phi_0 [graus]';
    vxlim = [0 175];
    FE = FE_ruido;
    KE = KE_ruido;
    FE1 = FE1_ruido;
    FE2 = FE2_ruido;
    case 2
        % --- tau
    load(pasta+"estimadoresFR_salto_freq_tau_n.mat");
    var="tau"; VEC=tau_n;
    sxlabel = '\tau_n [amostras]';
    vxlim = [48 432];
    KE = KFE1;
    case 3
        % ---- f1
    load(pasta+"estimadoresFR_salto_freq_F1.mat");
    VEC = F1_vec;
    var = "F1";
    sxlabel = 'f_1 [Hz]';
    vxlim = [58 62];
    KE = KFE;
    case 4
        % ---- altura do salto hf
        load(pasta+"estimadoresFR_salto_freq_kf.mat");
        VEC = kf_vec;
        sxlabel = 'h_f [Hz]';
        var = "hf";
        vxlim = [-3 3];
        KE = KFE;
        
    case 5
        % ---- frequencia de amostragem fs
        load(pasta+"estimadoresFR_salto_freq_Fs.mat");
        VEC = Fs_vec;
        sxlabel = 'f_s [Hz]';
        var = "fs";
        vxlim = [1000 20000];
        KE = KFE;
    case 6 
        % ---- tamanho da janela T
        load(pasta+"estimadoresFR_salto_freq_T.mat");
        VEC = Ncycles_vec;
        sxlabel = 'T [ciclos]';
        var = "T";
        vxlim = [3 30];
        KE = KFE;
     case 7 
        % ---- erro de estimação de tau
        load(pasta+"estimadoresFR_salto_freq_tau_error.mat");
        VEC = tau_n_error;
        sxlabel = '\epsilon [\Delta t]';
        var = "tau_error";
        vxlim = [-10 10];
        KE = KFE;
end

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
end

message = plot_type + " graficos criados na pasta "+ pasta;