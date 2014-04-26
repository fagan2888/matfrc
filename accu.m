
n = 1000;
my_error = zeros(1, 1000);
sys_error = zeros(1, 1000);

for i = 1:n
    test = randn(5);
    [zp, zl, zu] = zplu(test);
    [l, u] = lu(test);

    my_error(i) = norm(zl*zu - zp*test, 'fro');
    sys_error(i) = norm(l*u - test, 'fro');
end

fprintf('mean of my lu     : %f\n', mean(my_error));
fprintf('variance of my lu : %f\n', var(my_error));

fprintf('mean of matlab lu     : %f\n', mean(sys_error));
fprintf('variance of matlab lu : %f\n', var(sys_error));