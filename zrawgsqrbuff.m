
n = 2000;
m = 10;

qmyerrors = zeros(n, 1);
qstderrors = zeros(n, 1);

amyerrors = zeros(n, 1);
astderrors = zeros(n, 1);
for i=1:n
    
    test = randn(m, m);
    
    [zq, zr] = zrawgsqr(test);
    [q, r] = qr(test);
    
    %check the orthogonality of matrix zq and q
    qmyerrors(i) = norm(zq*(zq') - eye(m), 'fro');
    qstderrors(i) = norm(q*(q') - eye(m), 'fro');
    
    amyerrors(i) = norm(zq*zr - test, 'fro');
    astderrors(i) = norm(q*r - test, 'fro');
end

fprintf('Comparion of orthogonality of matrix Q\n')
fprintf('mean of my     norm(q*t(q) - eye(n) : %g\n', mean(qmyerrors));
fprintf('mean of matlab norm(q*t(q) - eye(n) : %g\n', mean(qstderrors));
fprintf('variance of my     norm(q*t(q) - eye(n) : %g\n', var(qmyerrors));
fprintf('variance of matlab norm(q*t(q) - eye(n) : %g\n', var(qstderrors));

fprintf('Comparison of QR decomposition\n');
fprintf('mean of my     norm(q*r - test) : %g\n', mean(amyerrors));
fprintf('mean of matlab norm(q*r - test) : %g\n', mean(astderrors));
fprintf('variance of my     norm(q*r - test) : %g\n', var(amyerrors));
fprintf('variance of matlab norm(q*r - test) : %g\n', var(astderrors));