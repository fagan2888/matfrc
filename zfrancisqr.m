function [H, U] = zfrancisqr(A)
%compute one of the step by implicitly shifted QR step

[m, n] = size(A);

if m ~= n
    error('support square matrix only')
end

m = n-1;

s = A(m, m) + A(n, n);
t = A(m, m)*A(n, n) - A(m, n)*A(n, m);

x = A(1,1)*A(1,1) + A(1,2)*A(2,1) - s*A(1,1) + t;
y = A(2,1)*(A(1,1) + A(2,2) - s);
z = A(2,1)*A(3,2);

for k=0:n-3
    
    [v, beta] = zhouse([x y z]');
    
    q = max([1 k]);
    
    %orthogonal transformation
    ot = (eye(3) - beta*v*(v'));
    
    A(k+1:k+3,q:n) = ot*A(k+1:k+3, q:n);
    
    r = min([k+4 n]);
    A(1:r, k+1:k+3) = A(1:r, k+1:k+3)*(ot');
    
    x = A(k+2, k+1);
    y = A(k+3, k+1);
    if k < n-3
        z = A(k+4, k+1);
    end
end

[v, beta] = zhouse([x y]');
ot = eye(2) - beta*v*(v');
A(n-1:n, n-2:n) = ot*A(n-1:n, n-2:n);
A(1:n, n-1:n) = A(1:n, n-1:n)*(ot');

H = A;
U = eye(n);