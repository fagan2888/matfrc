function [L, U] = zlu(A)
% ZLU - LU decomposition for matrix A
% work as gauss elimination

[m, n] = size(A);
if m ~= n 
    error('Error, current time only support square matrix');
end

L = zeros(n);
U = zeros(n);

for k = 1:n-1
    gauss_vector = A(:,k);
    gauss_vector(k+1:end) = gauss_vector(k+1:end) ./ gauss_vector(k);
    gauss_vector(1:k) = zeros(k,1);
    L(:,k) = gauss_vector;
    L(k,k) = 1;
    for l=k+1:n
        A(l,:) = A(l,:) - gauss_vector(l)*A(k,:);
    end
end    

U = A;