function [G] = zgaxpychol(A)
%cholesky decomposition for symmetric positive definite matrix
%the only requirement is matrix A: symmetric positive definite

[m, n] = size(A);

if m ~= n
    error('support square matrix only')
end

G = eye(n);

for k=1:n
    
    v = A(:,k);
    
    if k > 1
        v(:) = v(:) - G(:,1:k-1)*G(k,1:k-1)';
    end
    
    G(k:end, k) = v(k:end) / sqrt(v(k));
end

    