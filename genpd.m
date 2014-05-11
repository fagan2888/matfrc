function [A] = genpd(n)
%generate nXn positive definite matrix

temp = gensys(n);
[v, d] = eigs(-temp);
min_eigvalue = max(diag(d));

if min_eigvalue < 0
    A = temp;
else
    A = temp + (min_eigvalue + 0.5)*eye(size(temp));
end
