function [NW, NY] = zhousematrix(W, Y, v, beta)
%the product of n householder matrix can be formulated as Q = I + W*Y; for
%a new household matrix Q can be updated to Q';

[wm, wn] = size(W);
[ym, yn] = size(Y);

if wm ~= ym || wn ~= yn
    error('W and Y must have the same size')
end

[m, n] = size(v);
if m ~= wm
    error('the v and W, Y must have same dimension')
end

z = -beta*(eye(n) + W*(Y'))*v;
NW = [W z];
NY = [Y v];
