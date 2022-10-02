load('linconstr63.mat');
x0 = zeros(n,1);
for i = [7 16 17 20 22 23 26 30 40 44]
    x0(i) = 1;
end
A*x0-b