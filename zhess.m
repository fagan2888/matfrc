function [H, U] = zhess(A)
%for any matrix A, turn it into a upper hessenberg matrix by orthogonal 
%transformation

[m, n] = size(A);

if m ~= n
    error('support square matrix only')
end

H = A;
U = eye(n);

for k=1:n-2
    
    %compute the householder matrix
    [v, beta] = zhouse(H(k+1:end, k));
    temp_U = eye(n);
    
    temp_U(k+1:n,k+1:n) = eye(n-k) - beta*v*(v');
    
    H = temp_U*H;
    U = U * temp_U;
    
    %fprintf('after %d iteration\n', k);
    %disp(H);
end