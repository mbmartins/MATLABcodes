function [x,p,cost] = patv_Lp(y, d, lambda, p, E, Nit, mu0, mu1)
% [x,p,cost] = patv_Lp(y, d, lambda, p, E, Nit, mu0, mu1)
% Enhanced PATV: Simultaneous polynomial approximation and total variation
% filtering
%   Regularization : lambda * sum((abs(diff(x)) + E).^p);
%
% INPUT
%   y - noisy data
%   d - order of polynomial
%   lambda - regularization parameter
%   p, E - Lp norm
%   Nit - number of iterations
%   mu0, mu1 - Augmented Lagrangian parameters
%
% OUTPUT
%   x - TV component
%   p - polynomial component
%   cost - cost function history

% Ivan Selesnick
% Polytechnic Institute of New York University
% December 2011
% Reference: Polynomial Smoothing of Time Series with Additive Step Discontinuities
% I. W. Selesnick, S. Arnold, and V. R. Dantham

y = y(:);                   % convert to column vector
N = length(y);

n = (0:N-1)';
G = zeros(N,d);
for k = 1:d, G(:,k) = n.^k;  end        % exclude dc term (included in TV component)
G = orth(G);                            % orthogonalize cols of G (so that G' G = I)

e = ones(N-1,1);
D = spdiags([-e e],[0 1],N-1,N);        % D (sparse matrix)
A = mu0*(D'*D) + mu1*speye(N);          % mu0 D'D + mu1 I (sparse matrix)
D = @(x) diff(x);                       % D (operator)
DT = @(x) [-x(1); -diff(x); x(end)];    % D' (operator)
H = @(x) x - G * (G'*x);

x = zeros(N,1);                         % initializations
d0 = zeros(N-1,1);
d1 = zeros(N,1);

M = 15;					% M : number of outer iterations
cost = zeros(M,Nit);    % cost function history
b = lambda;             % initialize

for i = 1:M
    for k = 1:Nit
        u0 = soft(D(x)+d0, 0.5*b/mu0);
        tmp = x + d1;
        u1 = ((y + mu1*tmp) + G*(G'*(tmp-y))) / (1+mu1);
        x = A \ (mu0*DT(u0-d0) + mu1*(u1-d1));         % sparse system solve
        d0 = d0 - (u0-D(x));
        d1 = d1 - (u1-x);
        cost(i,k) = sum(abs(b .* D(x))) + sum(abs(H(x-y)).^2);
    end
    b = lambda * p * (abs(diff(x)) + E).^(p-1);
end
p = G * (G' * (y - x));
