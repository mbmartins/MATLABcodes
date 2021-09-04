function [V_pmu,freq] = PMU_FFT(PMU_sig, Fs)

        V_pmu_mag = 2*abs(fft(PMU_sig))/(length(PMU_sig));
        V_pmu_ang = angle(fft(PMU_sig));
        [vmax,imax] = max(V_pmu_mag);
        V_pmu = (V_pmu_mag(imax)/sqrt(2))*exp(i*V_pmu_ang(imax));
        bin = Fs/length(PMU_sig);
        freq = (imax-1)*bin;