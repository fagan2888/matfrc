function [Q, R] = zhouseqr(A)
%compute the QR decomposition of square matrix A
% A = Q*R; Q is orthogonal matrix and R is upper trigular matrix

[m, n] = size(A);

if m ~= n
    error('support squre matrix only')
end

temp = A;
Q = eye(n);

for k=1:n-1
    
    [v, beta] = zhouse(temp(k:end, k));
    reflector = eye(m-k+1) - beta*v*(v');
    temp(k:m, k:m) = reflector*temp(k:m, k:m);
    
    %construct orthogonal matrix
    t = eye(n);
    t(k:m, k:m) = reflector;
    Q = Q*(t');
    %
end

R = temp;