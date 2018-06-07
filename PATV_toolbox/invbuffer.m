function x = invbuffer(X, P, w)
% x = invbuffer(X, P)
% inverts the buffer function
%
% x = invbuffer(X, P, w)
% Multiplies each frame by window w
%
% Each column of X is one block of x
% P : overlap

% Ivan Selesnick
% Polytechnic Institute of New York University
% December 2011

[L, M] = size(X);
% L : length of block
% M : number of blocks

x = zeros(L+(M-1)*(L-P),1);

if nargin < 3    
    for i = 1:M
        x((i-1)*(L-P)+(1:L)) = x((i-1)*(L-P)+(1:L)) + X(:,i);
    end    
else    
    s = x; % zeros
    w = w(:);
    for i = 1:M
        x((i-1)*(L-P)+(1:L)) = x((i-1)*(L-P)+(1:L)) + w .* X(:,i);
        s((i-1)*(L-P)+(1:L)) = s((i-1)*(L-P)+(1:L)) + w;
    end
    x = x./s;    
end