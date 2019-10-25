function [x, p, cost, u, v,a] = patv_MM_modified(y, d, lambda, Nit)
% [x, p, cost] = patv_MM(y, d, lambda, Nit)
% PATV: Simultaneous polynomial approximation and total variation denoising.
% Algorithm based on MM algorithm.
%
% INPUT
%   y - noisy data
%   d - order of polynomial
%   lambda - regularization parameter
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

S = @(x) [0; cumsum(x)];                % S : cumulative sum 
rev = @(x) x(end:-1:1,:);               % rev : reverse elements of vector
ST = @(x) rev(cumsum(rev(x(2:end,:)))); % ST : S'
H = @(x) x - G * (G'*x);                % H : I - G G'
B = ST(G);
C = G'*G;

b = ST(H(y));                           % b : S'*H(y);
u = ones(N-1, 1);                       % Initialization
x = S(u);
Dx = u;                                 % Dx = D*x = D*S*u = I*u = u

for k = 1:Nit        
    Lam = spdiags(abs(Dx)/lambda, 0, N-1, N-1);
    R = STSinv + Lam;                           % size(R) : (N-1) x (N-1)
    invA = @(x) Lam * (x - R \ (Lam * x));      % invA : (N-1) -> (N-1)
    Q = C - B' * invA(B);                       % Q : size(Q) = (d+1) x (d+1)
    u = invA(b + B * (Q \ (B' * invA(b))));    
    x = S(u);
    Dx = u;                             % Dx = D*x = D*S*u = I*u = u
    cost(k) = lambda * sum(abs(Dx)) + 0.5 * sum(abs(H(x - y)).^2);
end

a = G'*(y - x); %obtida na base G
p = G * (G' * (y - x));
v = ST(H(y - x));
