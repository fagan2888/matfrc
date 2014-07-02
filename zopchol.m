function [G] = zopchol(A)
%cholesky decomposition based on rank-1 matrix update
%this problem can be formalized as A = [a b; b' B];
%then A = [ c d; d' D ]; and c = sqrt(a); d = b / c; and D = B - d*d' / a ;

[m, n] = size(A);
if m ~= n
    error('support square matrix only')
end

G = zeros(n);

for k=1:n
    
    G(k,k) = sqrt(A(k,k));
    G(k+1:end, k) = A(k+1:end, k) / G(k,k);
    
    %update matrix A
    for j = (k+1):n
        A(k+1:end,j) = A(k+1:end,j) - G(j,k)*G(k+1:end,k);
    end
end

    
