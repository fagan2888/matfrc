
n = 1500;
my_error = zeros(1, n);
sys_error = zeros(1, n);

for i = 1:n
    test = gensys(5);
    [zd, zl] = zldl(test);
    [l, d] = ldl(test);

    my_error(i) = norm(zl*zd*(zl') - test, 'fro');
    sys_error(i) = norm(l*d*(l') - test, 'fro');
end

fprintf('mean of my lu     : %g\n', mean(my_error));
fprintf('variance of my lu : %g\n', var(my_error));

fprintf('mean of matlab lu     : %g\n', mean(sys_error));
fprintf('variance of matlab lu : %g\n', var(sys_error));