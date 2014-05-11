function [A] = gensys(n)
% generate sysmetric matrix

A = randn(n);

for i=1:n
    for j=1:i
        A(j,i) = A(i,j);
    end
end