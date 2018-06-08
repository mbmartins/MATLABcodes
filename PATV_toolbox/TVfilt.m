function [x, cost] = TVfilt(y, lambda, Nit, mu1, mu2)
% [x, cost] = TVfilt(y, lambda, Nit, mu1, mu2)
% Total variation filtering (denoising) using ADMM
%
% INPUT
%   y - noisy signal
%   lambda - regularization parameter
%   Nit - number of iterations
%   mu1, mu2 - augmented Lagrangian parameters
%
% OUTPUT
%   x - filtered signal
%   cost - objective function history

% Ivan Selesnick
% Polytechnic Institute of New York University


y = y(:);                               % convert to column vector
cost = zeros(1, Nit);                   % objective function
N = length(y);
T = 0.5*lambda/mu2;

e = ones(N-1, 1);
Dmtx = spdiags([e -e], [0 1], N-1, N);  % D : first-order diff matrix (sparse)
F = mu1*speye(N) + mu2*(Dmtx'*Dmtx);    % F = mu1 I + mu2 D'*D (sparse)

D = @(x) diff(x);                       % D
DT = @(x) [-x(1); -diff(x); x(end)];    % D'

% initializations
d1 = zeros(N, 1);
d2 = zeros(N-1, 1);
x = zeros(N, 1);                        

for k = 1:Nit
    v1 = (y - d1 + mu1*x) / (1 + mu1);
    v2 = soft(D(x) + d2, T) - d2;    
    x = F \ (mu1*v1 + mu2*DT(v2));      % sparse system solver
    d1 = x - v1;
    d2 = D(x) - v2;
    cost(k) = sum(abs(x-y).^2) + lambda * sum(abs(D(x)));    
end

