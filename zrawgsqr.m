function [Q, R] = zrawgsqr(A)
%compute the QR decomposition based on raw gram-schmit approach

[m, n] = size(A);

if m ~= n
    error('support square matrix only');
end

Q = zeros(n, n);
R = zeros(n, n);

for k=1:n
    
    for j=1:k-1
        R(j, k) = Q(:,j)' * A(:, k);
        A(:, k) = A(:, k) - R(j, k)*Q(:,j);
    end
    
    R(k, k) = norm(A(:,k), 2);
    Q(:, k) = A(:, k) ./ R(k, k);
end
