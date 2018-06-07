function [x, p, p_blocks, cost, constr, r] = clopatv(y, L, P, deg, ES, Nit, mu0, mu1, lambda)
% [x, p, p_blocks, cost, constr] = clopatv(y, L, P, deg, ES, Nit, mu0, mu1)
% Simultaneous local polynomial approximation and total variation filtering
% - constrained formulation
%
% INPUT
%   y - noisy data
%   L - block length
%   P - overlapping (number of samples common to adjacent blocks)
%   deg - polynomial degree
%   ES - estimated noise sigma (standard deviation)
%   Nit - number of iterations
%   mu0, mu1 - augmented Lagrangian parameters
%
% OUTPUT
%   x - TV component
%   p - polynomial component
%   cost - cost function history
%   constr - constraint function history

% Use clopatv(..., lambda) for weighted TV minimization

y = y(:);                   % convert to column vector
cost = zeros(1,Nit);        % cost function history
constr = zeros(1,Nit);      % constraint function history
N = length(y);

if nargin < 9, lambda = 1; else lambda = lambda(:); end

G = zeros(L,deg+1);
for k = 0:deg, G(:,k+1) = (0:L-1)'.^k; end
G = orth(G);                                % orthogonalize cols of G

H = @(x) x - G * (G'*x);
buff = @(x) buffer(x, L, P, 'nodelay');

B = buff(1:N);
M = size(B,2);                              % M - number of blocks
s = invbuffer(buff(ones(1,N)), P);

e = ones(N-1,1);
D = spdiags([-e e],[0 1],N-1,N);            % D (sparse matrix)
A = mu0*(D'*D) + mu1*spdiags(s,0,N,N);      % mu0 D'D + mu S (sparse matrix)
D = @(x) diff(x);                           % D (operator)
DT = @(x) [-x(1); -diff(x); x(end)];        % D' (operator)

x = zeros(size(y));                         % initialization
xb = buff(x);
d0 = zeros(N-1,1);
d = zeros(size(xb));
Hyb = H(buff(y));

G1 = G(:,2:end);                            % omit dc term for first block
r = ES*sqrt(L*M);

G_ = zeros(N, M*(deg+1)-1);
for i = 1:M
     if i > 1
       G_((1:L)+(i-1)*(L-P) , (1:deg+1)+(i-1)*(deg+1)-1) = G;
    else
       G_((1:L)+(i-1)*(L-P) , 1:deg) = G1;
    end
end

F = (1/mu1)*eye((deg+1)*M-1) - (G_' * (A \ G_));    % F (sparse matrix solve)
G2_ = F \ G_';

for it = 1:Nit
    u0 = soft(D(x)+d0, 0.5*lambda/mu0);
    u = projball(H(xb) + d, Hyb, r);

    b = A \ (mu0*DT(u0-d0) + mu1*invbuffer(H(u-d), P));  % sparse matrix solve
    x = b + A \ ( G_ * (G2_ * b));                  % sparse matrix solve (unfortunate that G_ and G2_ are large...)
    xb = buff(x);
    
    d0 = d0 - (u0-D(x));
    d = d - (u-H(xb));
    
    cost(it) = sum(lambda.*abs(D(x)));
    constr(it) = sqrt(sum(sum(abs(H(buff(y-x))).^2)));
end
p_blocks = G * (G' * buff(y - x));

if 0
    % verify
    t = (Hyb - u).^2;
    sqrt(sum(t(:)))
    
    t = (H(buff(x)) - Hyb).^2;
    sqrt(sum(t(:)))
end

% compute a polynomial estimate
w = hamming(L);                     % column vector
p = invbuffer(p_blocks, P, w);
