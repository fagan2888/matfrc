function [P, Q, L, U] = zflu(A)
%full pivoted LU decomposition
%
% full pivoted LU decomposition

[m, n] = size(A);

if m ~= n
    error('current only support square matrix')
end

P = eye(n);
Q = eye(n);

for k=1:n-1
    
    %find the larget element in A(k:n,k:n)
    [max_value, row_index] = max(A(k:n, k:n));
    [max_value, col_index] = max(max_value);
    
    real_row = k-1 + row_index(col_index);
    real_col = k-1 + col_index;
    
    %exchange the row and column of matrix A
    
    if real_row ~= k
        A([k real_row],:) = A([real_row k], :);
        P([k real_row],:) = P([real_row k], :);
    end
    
    if real_col ~= k
        A(:, [k real_col]) = A(:, [real_col k]);
        Q(:, [k real_col]) = Q(:, [real_col k]);
    end
    
    if A(k, k) ~= 0
        rows = k+1:n;
        A(rows, k) = A(rows, k) ./ A(k, k);
        A(rows, rows) = A(rows, rows) - A(rows, k)*A(k, rows);
    end
end

L = tril(A);
for k=1:n
    L(k, k) = 1;
end
U = triu(A);

