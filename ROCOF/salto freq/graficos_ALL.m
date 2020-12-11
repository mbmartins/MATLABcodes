%geracao dos graficos
clc; clear all; c = ['k','b','c','r','m','g'];

% -----  influencia de phi0
load('salto_freq_phi0.mat');close all;
% ------------- ALL com ru�do de 60 dB -------------
for k = 1:6
    figALLphi0(k) = figure('Units','normalized','Position',[0 0 0.5 0.13]);
    subplot(1,3,1)
    plot(phi_n,zeros(1,length(phi_n)),'k--'); hold on;
    pleg(1) = plot(phi_n,FE_ruido(:,k),c(k)+"-"); hold on;
    shade(phi_n,FE_ruido(:,k) - FE_std(:,k),[c(k),'--'],phi_n,FE_ruido(:,k) + FE_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel('\phi_0 [graus]'); ylabel('FE [Hz]')
    xlim([0,175]); ylim([-0.06 0.06])
    grid on;
    subplot(1,3,2)
    plot(phi_n,zeros(1,length(phi_n)),'k--'); hold on;
    pleg(2) = plot(phi_n,KE_ruido(:,k),c(k)+"-"); hold on;
    shade(phi_n,KE_ruido(:,k) - kf_std(:,k),[c(k),'-.'],phi_n,KE_ruido(:,k) + kf_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
    xlabel('\phi_0 [graus]'); ylabel('K_fE [Hz]')
    xlim([0,175]); ylim([-1.1 0.2])
    grid on;
    subplot(1,3,3)
    plot(phi_n,zeros(1,length(phi_n)),'k--'); hold on;    
    pleg(3) = plot(phi_n,FE1_ruido(:,k),c(k)+"-"); hold on;
    shade(phi_n,FE1_ruido(:,k) - f1_std(:,k),[c(k),'-.'],phi_n,FE1_ruido (:,k) + f1_std(:,k),[c(k),'-.'], 'FillType', [1,2; 2,1]);
    xlabel('\phi_0 [graus]'); ylabel('FE_1 [Hz]')
    grid on;
    xlim([0,175]); ylim([-0.1 0.5])

savefig(figALLphi0(k),"salto_freq_ALL_phi0"+"_EF"+k)
saveas(figALLphi0(k),"salto_freq_ALL_phi0"+"_EF"+k+".png")
end
% --------------------------------------------------

% -----  influencia de tau_n
clc; clear all; c = ['k','b','c','r','m','g'];
load('estimadoresFR_salto_freq_tau_n.mat');close all;
VEC = tau_n; KFE = KFE1;
xlabels = '\tau_n [amostras]';
var = "tau";
% ------------- ALL com ru�do de 60 dB -------------
for k = 1:6
    figALLf1(k) = figure('Units','normalized','Position',[0 0 0.5 0.13]);
    subplot(1,3,1)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
    pleg(1) = plot(VEC,FE(:,k),c(k)+"-"); hold on;
    shade(VEC,FE(:,k) - FE_std(:,k),[c(k),'-.'],VEC,FE(:,k) + FE_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel(xlabels); ylabel('FE [Hz]')
    xlim([48,432]); ylim([-0.5 0.2])
    grid on;
    subplot(1,3,2)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
    pleg(2) = plot(VEC,KFE(:,k),c(k)+"-"); hold on;
    shade(VEC,KFE(:,k) - kf_std(:,k),[c(k),'--'],VEC,KFE(:,k) + kf_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel(xlabels); ylabel('K_fE [Hz]')
    xlim([48,432]); ylim([-1.5 5.])
    grid on;
    subplot(1,3,3)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;    
    pleg(3) = plot(VEC,FE1(:,k),c(k)+"-"); hold on;
    shade(VEC,FE1(:,k) - f1_std(:,k),[c(k),'--'],VEC,FE1(:,k) + f1_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel(xlabels); ylabel('FE_1 [Hz]')
    grid on;
    xlim([48,432]); ylim([-2.5 1.])
savefig(figALLf1(k),"salto_freq_ALL_"+var+"_EF"+k)
saveas(figALLf1(k),"salto_freq_ALL_"+var+"_EF"+k+".png")
end
% --------------------------------------------------

% influencia de F1
clc; clear all; c = ['k','b','c','r','m','g'];
load('estimadoresFR_salto_freq_F1.mat');close all;
VEC = F1_vec;
% ------------- ALL com ru�do de 60 dB -------------
for k = 1:6
    figALLf1(k) = figure('Units','normalized','Position',[0 0 0.5 0.13]);
    subplot(1,3,1)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
    pleg(1) = plot(VEC,FE(:,k),c(k)+"-"); hold on;
    shade(VEC,FE(:,k) - FE_std(:,k),[c(k),'-.'],VEC,FE(:,k) + FE_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel('f_1 [Hz]'); ylabel('FE [Hz]')
    ylim([-0.3 0.1])
    grid on;
    subplot(1,3,2)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
    pleg(2) = plot(VEC,KFE(:,k),c(k)+"-"); hold on;
    shade(VEC,KFE(:,k) - kf_std(:,k),[c(k),'--'],VEC,KFE(:,k) + kf_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel('f_1 [Hz]'); ylabel('K_fE [Hz]')
    ylim([-1.1 0.1])
    grid on;
    subplot(1,3,3)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;    
    pleg(3) = plot(VEC,FE1(:,k),c(k)+"-"); hold on;
    shade(VEC,FE1(:,k) - f1_std(:,k),[c(k),'--'],VEC,FE1(:,k) + f1_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel('f_1 [Hz]'); ylabel('FE_1 [Hz]')
    grid on;
    ylim([-0.3 0.6])
savefig(figALLf1(k),"salto_freq_ALL_F1"+"_EF"+k)
saveas(figALLf1(k),"salto_freq_ALL_F1"+"_EF"+k+".png")
end
% --------------------------------------------------

% -----  influencia de kf
clc; clear all; c = ['k','b','c','r','m','g'];
load('estimadoresFR_salto_freq_kf.mat');close all;
VEC = kf_vec;
xlabels = 'f_1 [Hz]';
var = "kf";
% ------------- ALL com ru�do de 60 dB -------------
for k = 1:6
    figALLf1(k) = figure('Units','normalized','Position',[0 0 0.5 0.13]);
    subplot(1,3,1)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
    pleg(1) = plot(VEC,FE(:,k),c(k)+"-"); hold on;
    shade(VEC,FE(:,k) - FE_std(:,k),[c(k),'-.'],VEC,FE(:,k) + FE_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel(xlabels); ylabel('FE [Hz]')
    xlim([48,432]); ylim([-0.3 0.1])
    grid on;
    subplot(1,3,2)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
    pleg(2) = plot(VEC,KFE(:,k),c(k)+"-"); hold on;
    shade(VEC,KFE(:,k) - kf_std(:,k),[c(k),'--'],VEC,KFE(:,k) + kf_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel(xlabels); ylabel('K_fE [Hz]')
    xlim([48,432]); ylim([-0.3 0.1])
    grid on;
    subplot(1,3,3)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;    
    pleg(3) = plot(VEC,FE1(:,k),c(k)+"-"); hold on;
    shade(VEC,FE1(:,k) - f1_std(:,k),[c(k),'--'],VEC,FE1(:,k) + f1_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel(xlabels); ylabel('FE_1 [Hz]')
    grid on;
    xlim([48,432]); ylim([-0.3 0.1])
savefig(figALLf1(k),"salto_freq_ALL_"+var+"_EF"+k)
saveas(figALLf1(k),"salto_freq_ALL_"+var+"_EF"+k+".png")
end
% --------------------------------------------------

% -----  influencia de fs
clc; clear all; c = ['k','b','c','r','m','g'];
load('estimadoresFR_salto_freq_Fs.mat');close all;
VEC = Fs_vec;
xlabels = 'f_s [Hz]';
var = "fs";
% ------------- ALL com ru�do de 60 dB -------------
for k = 1:6
    figALLf1(k) = figure('Units','normalized','Position',[0 0 0.5 0.13]);
    subplot(1,3,1)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
    pleg(1) = plot(VEC,FE(:,k),c(k)+"-"); hold on;
    shade(VEC,FE(:,k) - FE_std(:,k),[c(k),'-.'],VEC,FE(:,k) + FE_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel(xlabels); ylabel('FE [Hz]')
    xlim = [min(Fs_vec) max(Fs_vec)]; ylim([-0.3 0.1])
    grid on;
    subplot(1,3,2)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
    pleg(2) = plot(VEC,KFE(:,k),c(k)+"-"); hold on;
    shade(VEC,KFE(:,k) - kf_std(:,k),[c(k),'--'],VEC,KFE(:,k) + kf_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel(xlabels); ylabel('K_fE [Hz]')
    xlim = [min(Fs_vec) max(Fs_vec)];ylim([-1.1 0.1])
    grid on;
    subplot(1,3,3)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;    
    pleg(3) = plot(VEC,FE1(:,k),c(k)+"-"); hold on;
    shade(VEC,FE1(:,k) - f1_std(:,k),[c(k),'--'],VEC,FE1(:,k) + f1_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel(xlabels); ylabel('FE_1 [Hz]')
    grid on;
    xlim = [min(Fs_vec) max(Fs_vec)];ylim([-0.3 0.6])
savefig(figALLf1(k),"salto_freq_ALL_"+var+"_EF"+k)
saveas(figALLf1(k),"salto_freq_ALL_"+var+"_EF"+k+".png")
end
% --------------------------------------------------

% -----  influencia de T
clc; clear all; c = ['k','b','c','r','m','g'];
load('estimadoresFR_salto_freq_T.mat');close all;
VEC = Ncycles_vec;
xlabels = 'T [ciclos]';
var = "T";
% ------------- ALL com ru�do de 60 dB -------------
for k = 1:6
    figALLf1(k) = figure('Units','normalized','Position',[0 0 0.5 0.13]);
    subplot(1,3,1)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
    pleg(1) = plot(VEC,FE(:,k),c(k)+"-"); hold on;
    shade(VEC,FE(:,k) - FE_std(:,k),[c(k),'-.'],VEC,FE(:,k) + FE_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel(xlabels); ylabel('FE [Hz]')
    xlim = [min(Ncycles_vec) max(Ncycles_vec)]; ylim([-0.3 0.1])
    grid on;
    subplot(1,3,2)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;
    pleg(2) = plot(VEC,KFE(:,k),c(k)+"-"); hold on;
    shade(VEC,KFE(:,k) - kf_std(:,k),[c(k),'--'],VEC,KFE(:,k) + kf_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel(xlabels); ylabel('K_fE [Hz]')
    xlim = [min(Ncycles_vec) max(Ncycles_vec)];ylim([-1.1 0.1])
    grid on;
    subplot(1,3,3)
    plot(VEC,zeros(1,length(VEC)),'k--'); hold on;    
    pleg(3) = plot(VEC,FE1(:,k),c(k)+"-"); hold on;
    shade(VEC,FE1(:,k) - f1_std(:,k),[c(k),'--'],VEC,FE1(:,k) + f1_std(:,k),[c(k),'--'], 'FillType', [1,2; 2,1]);
    xlabel(xlabels); ylabel('FE_1 [Hz]')
    grid on;
    xlim = [min(Ncycles_vec) max(Ncycles_vec)];ylim([-0.3 0.6])
savefig(figALLf1(k),"salto_freq_ALL_"+var+"_EF"+k)
saveas(figALLf1(k),"salto_freq_ALL_"+var+"_EF"+k+".png")
end
% --------------------------------------------------