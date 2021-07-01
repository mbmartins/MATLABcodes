function pasta = gera_dados_MC(MCiter,h_f)
%gerar novamente os dados
%MCiter = 10000;
%h_f = 1; % salto de 1 Hz
pasta = "figuras_hf_"+h_f+"_"+MCiter+"\";
%phi_0
file_phi0 = teste_MedSF_phi0(MCiter,h_f,pasta);
%tau
file_tau = teste_MedSF_tau(MCiter,h_f,pasta);
%F1
file_F1 = teste_MedSF_F1(MCiter,h_f,pasta);
%hf
file_hf = teste_MedSF_hf(MCiter,h_f,pasta);
%fs
file_Fs = teste_MedSF_Fs(MCiter,h_f,pasta);
%T
file_T = teste_MedSF_T(MCiter,h_f,pasta);