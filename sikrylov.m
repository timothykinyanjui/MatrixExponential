function [f,decay] = sikrylov(A,b,m,shift)
% SIKRYLOV  Shift-and-invert Arnoldi method for expm(A)*b.
% Stefan Guettel, 2016.

if nargin < 4, shift = 1; end

B = A - shift*speye(size(A));
V = zeros(length(b),m); H = zeros(m);
beta = norm(b);
w = b/beta;
V(:,1) = w;
for c = 1:m,
    w =  B\V(:,c);
    
    % MGS
    for r = 1:c,
        H(r,c) = V(:,r)'*w;
        w = w - V(:,r)*H(r,c);
    end
    
    % CGS
    %H(1:c,c) = V(:,1:c)'*w;
    %w = w - V(:,1:c)*H(1:c,c);
    
    if c<m,
        H(c+1,c) = norm(w);
        V(:,c+1) = w/H(c+1,c);
    end
end
I = eye(m);
%E = expm(inv(H)); E1 = E(:,1);
E = expm(H\I); E1 = E(:,1); 
%[W,D] = eig(H); E1 = W*(exp(1./diag(D)).*(W\I(:,1))); % faster but unstable

decay = (beta*exp(shift))*E1;
f = V*decay;

