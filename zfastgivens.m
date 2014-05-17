function [alpha, beta, type, nd] = zfastgivens(x, d)
%compute parameter for fast givens rotation

%sanity check
[m, n] = size(x);
if m ~= 2
    error('size(x) must be [2, 1]')
end

if abs(x(2)) > 1e-5
    
    alpha = -1*x(1) ./ x(2);
    beta = -1*alpha*d(2) ./ d(1);
    gamma = -1*alpha*beta;
    
    if abs(gamma) <= 1
        type = 1;
        tau = d(1);
        d(1) = (1 + gamma)*d(2);
        d(2) = (1 + gamma)*tau;
    else
        type = 2;
        alpha = 1 ./ alpha;
        beta = 1 ./ beta;
        gamma = 1 ./ gamma;
        d(1) = (1 + gamma)*d(1);
        d(2) = (1 + gamma)*d(2);
    end
else
    type = 2;
    alpha = 0;
    beta = 0;
end

nd = d;