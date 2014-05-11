
n = 1500;
my_error = zeros(1, n);
sys_error = zeros(1, n);

for i = 1:n
    test = genpd(5);
    [zg] = zopchol(test);
    [l, u] = lu(test);

    my_error(i) = norm(zg*(zg') - test, 'fro');
    sys_error(i) = norm(l*u - test, 'fro');
end

fprintf('mean of my lu     : %g\n', mean(my_error));
fprintf('variance of my lu : %g\n', var(my_error));

fprintf('mean of matlab lu     : %g\n', mean(sys_error));
fprintf('variance of matlab lu : %g\n', var(sys_error));