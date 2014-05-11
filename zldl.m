function [D, L] = zldl(A)
%A = L*D*L' another version of LU decomposition for matrix A

[m, n] = size(A);

if m ~= n
    error('support square matrix only')
end

L = eye(n);
d = zeros(n,1);

for k=1:n
    v = zeros(n,1);
    
    for j=1:k-1
        v(j) = L(k, j)*d(j);
    end
    
    v(k) = A(k,k) - L(k, 1:k-1)*v(1:k-1);
    
    d(k) = v(k);
    
    L(k+1:end, k) = (A(k+1:end,k) - A(k+1:end, 1:k-1)*v(1:k-1)) / v(k);
end

D = diag(d);