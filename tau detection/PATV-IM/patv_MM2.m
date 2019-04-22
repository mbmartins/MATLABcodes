function [x, p, cost, u, v] = patv_MM(y, d, phi, wfun, Nit, u_init)
% [x, p, cost] = patv_MM2(y, d, phi, wfun, Nit)
% PATV: Simultaneous polynomial approximation and total variation
% denoising (implemented using MM algorithm)
%
% INPUT
%   y - noisy data
%   d - order of polynomial
%   phi - penalty function
%   wfun - function wfun(x) = x / phi'(x)
%   Nit - number of iterations
%
% OUTPUT
%   x - TV component
%   p - polynomial component
%   cost - cost function history

% Ivan Selesnick, selesi@poly.edu
% Polytechnic Institute of New York University
% November 2012
% Reference: 'Simultaneous least square polynomial approximation and total variation denoising'
% ICASSP 2013

y = y(:);                   % convert to column vector
cost = zeros(1,Nit);        % cost function history
N = length(y);

G = orth(bsxfun(@power, (0:N-1)', 0:d));

e = ones(N-1,1);
STSinv = spdiags([-e 2*e -e], [-1 0 1], N-1, N-1);
STSinv(1,1) = 1;

S = @(x) [0; cumsum(x)];            % S : cumulative sum
rev = @(x) x(end:-1:1,:);
ST = @(x) rev(cumsum(rev(x(2:end,:))));

H = @(x) x - G * (G'*x);
B = ST(G);
C = G'*G;
b = ST(H(y));                       % b : S' H y;

if nargin < 6
    u = ones(N-1, 1);               % initialization
else
    u = u_init;
end
Dx = u;                             % Dx = D*x = D*S*u = I*u = u

for k = 1:Nit
    Lam = spdiags(wfun(Dx), 0, N-1, N-1);    
    R = STSinv + Lam;                           % size(R) : (N-1) x (N-1)    
    invA = @(x) Lam * (x - R \ (Lam * x));      % invA : (N-1) -> (N-1)    
    Q = C - B' * invA(B);                       % Q : size(Q) = (d+1) x (d+1)    
    u = invA(b + B * (Q \ (B' * invA(b))));    
    x = S(u);
    Dx = u;                         % Dx = D*x = D*S*u = I*u = u    
    cost(k) = sum(phi(Dx)) + 0.5 * sum(abs(H(x - y)).^2);    
end

p = G * (G' * (y - x));
v = ST(H(y - x));