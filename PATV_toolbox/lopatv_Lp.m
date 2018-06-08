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
% OUTPUT
%   x - TV component
%   p - local polynomial component
%   cost - cost function history
%
% Number of blocks = (length(y)-L)/(L-P)+1
% If this is not an integer, then input signal y will be truncated.

% Ivan Selesnick
% Polytechnic Institute of New York University
% December 2011
% Reference: Polynomial Smoothing of Time Series with Additive Step Discontinuities
% I. W. Selesnick, S. Arnold, and V. R. Dantham

y = y(:);                               % convert to column vector
N = length(y);

M = (N-L)/(L-P);                        % M : number of blocks - 1
if M > floor(M)
    N = floor(M)*(L-P)+L;
    y = y(1:N);
    fprintf('Note in lopatv.m: The input signal will be truncated down to length %d\n',N)
    fprintf('so it is consistent with the block length %d and overlap %d.\n', L, P);
end

cost = zeros(1,Nit);                    % cost function history
G = zeros(L,deg+1);
for k = 0:deg, G(:,k+1) = (0:L-1)'.^k; end
G = orth(G);                            % orthogonalize columns of G (so that G' G = I)
H = @(x) x - G * (G'*x);
buff = @(x) buffer(x, L, P, 'nodelay');
s = invbuffer(buff(ones(1,N)), P);

e = ones(N-1,1);
D = spdiags([-e e],[0 1],N-1,N);                % D (sparse matrix)
A = mu0*(D'*D) + mu*spdiags(s,0,N,N);           % mu0 D'D + mu S (sparse matrix)
D = @(x) diff(x);                               % D (operator)
DT = @(x) [-x(1); -diff(x); x(end)];            % D' (operator)
x = zeros(size(y));                             % initialization
yb = buff(y);
xb = buff(x);
d0 = zeros(N-1,1);
d = zeros(size(xb));

M = 5;					% M : number of outer iterations (increase, if nec.)
cost = zeros(M,Nit);    % cost function history
b = lambda;             % initialize

for i = 1:M
    for k = 1:Nit
        u0 = soft(D(x)+d0, 0.5*b/mu0);
        tmp = xb + d;
        u = ((yb + mu*tmp) + G*(G'*(tmp-yb))) / (1+mu);
        x = A \ (mu0*DT(u0-d0) + mu*invbuffer(u-d, P));    % sparse system solve
        xb = buff(x);
        d0 = d0 - (u0-D(x));
        d = d - (u-xb);
        cost(i,k) = sum(abs(b .* D(x))) + sum(sum(abs(H(buff(y-x))).^2));
    end
    b = lambda * pow * (abs(diff(x)) + E).^(pow-1);
    % plot(1./b,'k'),  drawnow
end
p_blocks = G * (G' * buff(y - x));

w = hamming(L);                     % column vector
p = invbuffer(p_blocks, P, w);      % compute local polynomial estimate
