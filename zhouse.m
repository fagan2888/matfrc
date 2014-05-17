function [v, beta] = zhouse(x)
%generate householder reflector vector v, so that P=(I - 2*(v*v')/(v'*v)
% Px = epison*eye(:,1), good computation methods

[m, n] = size(x);

if n ~= 1
    error('only calculate householder reflector for vector not matrix')
end

delta = (x(2:end)')*x(2:end);

v = zeros(m, 1);
v(1) = 1;
v(2:end) = x(2:end);

if abs(delta) < 1e-8
    beta = 0;
else
    mu = sqrt( x(1)^2 + delta );
    
    if x(1) <= 0
        v(1) = x(1) - mu;
    else
        v(1) = -1*(delta ./ (x(1) + mu));
    end
    
    beta = 2 * (v(1)^2) ./ (delta + v(1)^2);
    v = v ./ v(1);
end