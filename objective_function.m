%% Aggregate score functions for different $k$, $k=1,3,5,8,10$

clear all
clc
close all

warning all

z = [45 8 22 91 15 82 53 7 44 99];

zsort = sort(z,'descend');

f={}

for k = [1 3 5 8 10]

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

tau1 = 1;
tau2 = 0.505;


alpha0 = 80;

alpha = @(t) alpha0/((t+1)^tau1); %fast_scale
beta = @(t) beta0/((t+1)^tau2);

p = (n-k)/n + 1/(2*n);

rho_p = @(x) (p-1)*x*(x<=0) + p*x*(x>0);

f{k} = @(x) rho_p(w(1)-x) + rho_p(w(2)-x) + rho_p(w(3)-x) + rho_p(w(4)-x) + rho_p(w(5)-x) + rho_p(w(6)-x) + rho_p(w(7)-x) + rho_p(w(8)-x) + rho_p(w(9)-x) + rho_p(w(10)-x);

end

%%
figure
box on
hold on
fplot(f{1},[0,110],'-s')
fplot(f{3},[0,110],'-x')
fplot(f{5},[0,110],'-^')
fplot(f{8},[0,110],'-*')
fplot(f{10},[0,110],'-o')

h=legend('$k=1$','$k=3$','$k=5$','$k=8$','$k=10$');
set(h,'Interpreter','latex');

xlabel('$\xi$','interpreter','latex')
ylabel('$f^k(\xi)$','interpreter','latex')
axis([0,110,0,700])