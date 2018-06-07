function [x, p, cost] = lopatv_Lp(y, L, P, deg, lambda, Nit, mu0, mu, pow, E)
% [x, p, cost] = lopatv_Lp(y, L, P, deg, lambda, Nit, mu0, mu, pow, E)
% LoPATV_Lp: Enhanced local polynomial approximation and total variation filtering
% (sliding window with overlapping) with Lp norm
%
% INPUT
%   y - noisy data
%   L - block length
%   P - overlapping (number of samples common to adjacent blocks)
%   deg - polynomial degree (1, 2, 3)
%   lambda - regularization parameter
%   Nit - number of iterations
%   mu - augmented Lagrangian parameters
%   pow - power (Lp norm)
%   E - small number
%
%
% OUTPUT
%   x - TV component
%   p - local polynomial component
%   cost - cost function history

% Ivan Selesnick
% Polytechnic Institute of New York University
% December 2011
% Reference: Polynomial Smoothing of Time Series with Additive Step Discontinuities
% I. W. Selesnick, S. Arnold, and V. R. Dantham

[x, p, cost] = lopatv(y, L, P, deg, lambda, Nit, mu0, mu);

for it = 1:5
    b = lambda * pow * (abs(diff(x)) + E).^(pow-1);
    % plot(1./b,'k'),  drawnow
    lambda_new = b;    
    [x, p, cost] = lopatv(y, L, P, deg, lambda_new, Nit, mu0, mu);    
end

