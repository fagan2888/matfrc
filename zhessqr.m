function [NH] = zhessqr(H)
% perform QR algorithm on upper hessenberg matrix
% firstly, we need to verity this is a hessberg matrix


[m, n] = size(H);

if m ~= n
    error('error, support square matrix only')
end

NH = H;

c = zeros(1, n-1);
s = zeros(1, n-1);

for k=1:n-1
    %compute gives rotation at first
    [c(k), s(k)] = zgivens(NH(k, k), NH(k+1, k));
    p = [c(k) s(k); -s(k) c(k)];
    NH(k:k+1, k:n) = (p')*NH(k:k+1, k:n);
    %fprintf('after %d iteration\n', k);
    %disp(NH);
end

for k=1:n-1
    p = [c(k) s(k); -s(k) c(k)];
    NH(1:k+1, k:k+1) = NH(1:k+1,k:k+1)*p;
end