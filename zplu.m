function [P, L, U] = zplu(A)
% pivoted LU decompositon P*A = L*U

[m, n] = size(A);

if m ~= n
    error('zplu:test', 'current time only support square matrix');
end

P = eye(n);
L = zeros(n, n);

for k = 1:n-1

    %find the largest element in k column of A from row k to n
    [max_value, max_index] = max(A(k:end, k));
    
    max_index = max_index + k - 1;
    if max_index ~= k
        A([k max_index], :) = A([max_index k], :);
        P([k max_index], :) = P([max_index k], :);
        L([k max_index], :) = L([max_index k], :);
    end
    
    if A(k,k) ~= 0
        gauss_vector = A(:,k);
        gauss_vector(k+1:end) = gauss_vector(k+1:end) ./ gauss_vector(k);
        gauss_vector(1:k) = zeros(k,1);
        L(:,k) = gauss_vector;
        L(k, k) = 1;
    
        for l=k+1:n
            A(l,:) = A(l,:) - gauss_vector(l)*A(k,:);
        end
    end
end
U = triu(A);
