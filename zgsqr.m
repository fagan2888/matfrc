function [Q, R] = zgsqr(A)
%compute QR decomposition based on intuitive way

[m, n] = size(A);

if m ~= n
    error('support square matrix only')
end

Q = zeros(n, n);
R = zeros(n, n);

for k = 1:n
    
    R(k, k) = norm(A(:,k), 2);
    Q(:, k) = A(:,k) ./ R(k, k);
    
    for j=k+1:n
        R(k, j) = Q(:,k)'*A(:,j);
        A(:, j) = A(:, j) - Q(:,k)*R(k,j);
    end
end