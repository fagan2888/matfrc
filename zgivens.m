function [c, s] = zgivens(a, b)
%compute givens rotation for vector [a b]'

if abs(b) < 1e-8
    c = 1; s = 0;
else
    if abs(b) > abs(a)
        tau = -1*a ./ b ;
        s = 1 ./ sqrt(1 + tau.^2);
        c = s * tau;
    else
        tau = -1*b ./ a;
        c = 1 ./ sqrt(1 + tau.^2);
        s = c * tau;
    end
end