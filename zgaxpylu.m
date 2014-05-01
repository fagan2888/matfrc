function [L, U] = zgaxpylu(A)

%calculate LU decomposition based on Gaxpy operation
%the same way as zlu.m but differnt approach

[m, n] = size(A);

if m ~= n
    error('current support square matrix only')
end

L = eye(n);
U = zeros(n);

for k = 1:n
    
    v = zeros(n, 1);
    if k == 1
        v(k:end) = A(k:end, k);
    else
        %solve triangular equation
        U(1:k-1,k) = L(1:k-1,1:k-1)\A(1:k-1,k);
        v(k:end) = A(k:end, k) - L(k:end, 1:k-1)*U(1:k-1,k);
    end
    
    if k < n
        L(k+1:end,k) = v(k+1:end)/v(k);
        L(k, k) = 1;
    end
    
    U(k, k) = v(k);
end
    
