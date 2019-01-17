function [f,A,ph] = ipfft(x,fs,p,window,model)

N = size(x,2);
T = 1/fs;

xw = x.*window;

Xfft = fft(xw);
Wfft = fft(window);
[Y,km] = max(Xfft(1:N/2));

alfa = abs(Xfft(km-1)); %eq (6)
beta = abs(Y);          %eq (7)
gamma = abs(Xfft(km+1));%eq (8)

%ip-dft classic
if strcmp(model,'parabola')
   K_ = km-1 + 0.5*(alfa-gamma)/(alfa-2*beta+gamma); %eq (11)
   X_ = beta - 0.125*(alfa-gamma)^2/(alfa-2*beta+gamma); %eq (12)
elseif strcmp(model,'log')
%LQIFFT
  falfa = log(alfa); fbeta = log(beta); fgamma = log(gamma);
  K_ = km-1 + 0.5*(falfa-fgamma)/(falfa-2*fbeta+fgamma); %eq (13)
  X_ = exp(fbeta - 0.125*(falfa-fgamma)^2/(falfa-2*fbeta+fgamma)); %eq (14)
elseif strcmp(model,'power')
%XQIFFT
    falfa = alfa^p; fbeta = beta^p; fgamma = gamma^p;
    K_ = km-1 + 0.5*(falfa-fgamma)/(falfa-2*fbeta+fgamma); %eq (13)
    X_ = (fbeta - 0.125*(falfa-fgamma)^2/(falfa-2*fbeta+fgamma))^(1/p); %eq (14)
end
%results
f = K_/(N*T);
A = 2*abs(X_)/Wfft(1);
ph = angle(Xfft(km-1));