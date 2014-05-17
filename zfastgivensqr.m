function [Q, R] = zfastgivensqr(A)
%QR decomposition based on fastgivens

[m, n] = size(A);

if m ~= n
    error('only support square matrix')
end

d = ones(m,1);
M = eye(n);

for k=1:n
    for i=m:-1:k+1
        
        [alpha, beta, type, nd] = zfastgivens(A([i-1 i], k), d([i-1 i]));
        d([i-1 i]) = nd;
        
        %construct M matrix explictly
        tempM = eye(n);
        
        if type == 1
            tempM(i-1, i-1) = beta;
            tempM(i-1, i) = 1;
            tempM(i, i-1) = 1;
            tempM(i, i) = alpha;
        else
            tempM(i-1, i-1) = 1;
            tempM(i-1, i) = alpha;
            tempM(i, i-1) = beta;
            tempM(i, i) = 1;
        end

        tx = A;
        A = (tempM')*A;
        
        %update dial matrix
        M = M*tempM;
        tempQ = (M*diag(1 ./ sqrt(d)));
    end
end

R = A;
Q = (M*diag(1 ./ sqrt(d)))';