function filename = Robustez_MedSF(phi_n,tau_vec,MCiter,h_f,filedisc)
%testes para robustez do MedSF aos saltos de fase e magnitude
SNR = 60;
%fixed parameters
F0 = 60.0; %nominal frequency
F1 = 60; %fundamental frequency
Fs =4800;
NCycles = 6;
T = NCycles/F0;
NSamples = floor(NCycles*Fs/F0);
nbits = 16;

% MC loop
% for k = 1:MCiter
%     k
%baseline
%k_a = 0; k_x = 0; k_f = 0;
%[FEraw_zero,fE1raw_zero,fE2raw_zero,kfEraw_zero, riraw_zero, draw_zero] = MC_estimation_Robustez(MCiter,F0,F1,Fs,phi_n(k),NCycles,tau_vec,SNR,k_a, k_x,k_f,nbits);
%salto fase
k_a = 10; k_x = 0; k_f = 0;
[FEraw_fase,fE1raw_fase,fE2raw_fase,kfEraw_fase, riraw_fase, draw_fase] = MC_estimation_Robustez(MCiter,F0,F1,Fs,phi_n,NCycles,tau_vec,SNR,k_a, k_x,k_f,nbits);
kfest_fase = kfEraw_fase + k_f;
%salto mag
k_a = 0; k_x = 0.1; k_f = 0;
[FEraw_mag,fE1raw_mag,fE2raw_mag,kfEraw_mag, riraw_mag, draw_mag] = MC_estimation_Robustez(MCiter,F0,F1,Fs,phi_n,NCycles,tau_vec,SNR,k_a, k_x,k_f,nbits);
kfest_mag = kfEraw_mag + k_f;
%salto freq
k_a = 0; k_x = 0; k_f = h_f;
[FEraw_freq,fE1raw_freq,fE2raw_freq,kfEraw_freq, riraw_freq, draw_freq] = MC_estimation_Robustez(MCiter,F0,F1,Fs,phi_n,NCycles,tau_vec,SNR,k_a, k_x,k_f,nbits);
kfest_freq = kfEraw_freq + k_f;
% end
filename = "Robustez\Robustez"+MCiter+ "_hf" +h_f*10 + filedisc;
save(filename)

