%% Communication graph with n=10 and empirical CDF with k=3

clear all
clc
close all

z = [45 8 22 91 15 82 53 7 44 99];

zsort = sort(z,'descend');

k = 3;

top_k_truth = (z'>=zsort(k));

n=length(z);

A = zeros(n,n);

A(1,[2 4]) = 1;
A(2,[1 5 3]) = 1;
A(3,[2 5]) = 1;
A(4,[8 6 1]) = 1;
A(5,[2 3 7 10]) = 1;
A(6,[4 8 9]) = 1;
A(7,[5 10]) = 1;
A(8,[4 6]) = 1;
A(9,[6 10]) = 1;
A(10,[7 5 9]) = 1;

D = diag(sum(A));

L = D - A;

lambda=eig(L);

beta0 = 2/(lambda(2)+lambda(n));

w = z';

f1 = @(x) (w(1)<=x);
f2 = @(x) (w(2)<=x);
f3 = @(x) (w(3)<=x);
f4 = @(x) (w(4)<=x);
f5 = @(x) (w(5)<=x);
f6 = @(x) (w(6)<=x);
f7 = @(x) (w(7)<=x);
f8 = @(x) (w(8)<=x);
f9 = @(x) (w(9)<=x);
f10 = @(x) (w(10)<=x);

F = @(x) (f1(x)+f2(x)+f3(x)+f4(x)+f5(x)+f6(x)+f7(x)+f8(x)+f9(x)+f10(x))/10;

%%
figure
plot(graph(A))

figure
fplot(F,[0,110])
hold on

p = (n-k)/n + 1/(2*n);
plot([0:1:110],p*ones(1,111),'--')
xlabel('$\xi$','interpreter','latex')
ylabel('$\widehat{F}(\xi ; \mathcal{D})$','interpreter','latex')