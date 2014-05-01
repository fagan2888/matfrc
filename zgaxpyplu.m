function [P, L, U] = zgaxpyplu(A)
%compute LU decomposition based on Gaxpy operation with pivoting
%aimed at improve the stability of LU decomposition

[m, n] = size(A);
if m ~= n
    error('current support square matrix only')
end

P = eye(n);
L = eye(n);
U = zeros(n);

for k = 1:n
    
    v = zeros(n, 1);
    
    if k == 1
        v(k:end) = A(k:end, 1);
    else
        U(1:k-1, k) = L(1:k-1, 1:k-1) \ A(1:k-1, k);
        v(k:n) = A(k:n, k) - L(k:n, 1:k-1)*U(1:k-1, k);
    end
    
    %find the largest element in v(k:end)
    [max_value, max_index] = max(v(k:end));
    max_index = max_index + k - 1;
    
    if max_index ~= k
        %exchange the max_index-th row and k-th row
        A([k max_index], k+1:n) = A([max_index k], k+1:n);
        %exchange the max_index-th row and k-th row in L
        if k > 1
            L([k max_index], 1:k-1) = L([max_index k], 1:k-1);
        end
        P([k max_index], :) = P([max_index k], :);
        v([k max_index]) = v([max_index k]);
    end
    
    if (abs(v(k)) > 1e-6) && (k < n)
        v(k+1:end) = v(k+1:end) / v(k);
        L(k+1:end, k) = v(k+1:end);
    end
    
    U(k, k) = v(k);
end