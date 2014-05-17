function [Q, R] = zgivensqr(A)
%compute QR decomposition based on givens rotation

[m, n] = size(A);

if m ~= n
    error('support square matrix')
end

Q = eye(n);
temp = A;

for k=1:n-1
    for j = m:-1:k+1
        [c, s] = zgivens(temp(j-1,k), temp(j,k));

        %update the original matrix
        for l=k:n
            tau1 = temp(j-1, l);
            tau2 = temp(j, l);
            temp(j-1, l) = c*tau1 - s*tau2;
            temp(j, l) = s*tau1 + c*tau2;
        end
        
        %update the orthogonal matrix
 
        for l=1:n
            tau1 = Q(j-1, l);
            tau2 = Q(j, l);
            Q(j-1, l) = c*tau1 - s*tau2;
            Q(j, l) = s*tau1 + c*tau2;
        end       

    end
end

R = temp;