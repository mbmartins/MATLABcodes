%teste randn

n = 1:10000;
a = (randn(1));
mean(a)

for n=1:10000
    b(n) = randn(1,1);
end
mean(b)

figure
histfit(b,20,'normal')