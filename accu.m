
n = 1500;
my_error = zeros(1, 1000);
sys_error = zeros(1, 1000);

for i = 1:n
    test = randn(5);
    [zp, zq, zl, zu] = zflu(test);
    [l, u] = lu(test);

    my_error(i) = norm(zl*zu - zp*test*zq, 'fro');
    sys_error(i) = norm(l*u - test, 'fro');
end

fprintf('mean of my lu     : %g\n', mean(my_error));
fprintf('variance of my lu : %g\n', var(my_error));

fprintf('mean of matlab lu     : %g\n', mean(sys_error));
fprintf('variance of matlab lu : %g\n', var(sys_error));