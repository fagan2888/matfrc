function [L, D, U] = zldu(A)
%LDU decomposition of square matrix A. The first step for Cholesky
%decomposition

[m, n] = size(A);
if m ~= n
    error('support square matrix only')
end

L = eye(n);
U = eye(n);
d = zeros(n,1);

for k=1:n
    
    v = zeros(n, 1);
    if k == 1
        v(k:end) = A(k:end, k);
    else
        m = L(1:k-1, 1:k-1) \ A(1:k-1, k);
        for j = 1:k-1
            U(j, k) = m(j) / d(j);
        end
        
        v(k:end) = A(k:end, k) - L(k:end, 1:k-1)*m(:);
    end
    
    d(k) = v(k);
    
    if k < n
        L(k+1:end, k) = v(k+1:end)/v(k);
    end
    
end

D = diag(d);