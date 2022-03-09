%teste 

tau_nA = 60;
NA = 120;
fsf1 = 80;

for m = 1:120
    tau_nB = tau_nA + m*fsf1
    NB(m) = round(tau_nB*(NA/tau_nA))
end