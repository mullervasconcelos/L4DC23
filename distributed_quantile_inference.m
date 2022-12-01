%% istributed sample quantile inference

clear all
clc
close all


z = [45 8 22 91 15 82 53 7 44 99];

zsort = sort(z,'descend');

gap=1;

J=[];


k=1;

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

tau1 = 1;
tau2 = 0.505;

alpha0 = 80;

alpha = @(t) alpha0/((t+1)^tau1); %fast_scale
beta = @(t) beta0/((t+1)^tau2);

p = (n-k)/n + 1/(2*n);

std = sqrt(1); %% standard deviation of communication noise

w = z';

t = 0;

while t<=1e5
    
    v = std*randn(n);
    
    v = sum(A.*v)';
    
    g = [(w(1)>=z(1)) - p; (w(2)>=z(2)) - p; (w(3)>=z(3)) - p; (w(4)>=z(4)) - p; (w(5)>=z(5)) - p; (w(6)>=z(6)) - p; (w(7)>=z(7)) - p; (w(8)>=z(8)) - p; (w(9)>=z(9)) - p; (w(10)>=z(10)) - p];
    
    wplus = (eye(n)-beta(t)*L)*(w-alpha(t)*g) + beta(t)*v;
     
    
    w = wplus;
    
    J = [J norm(w-zsort(k))/n];
    
    threshold = w - gap/2;
    
    t = t+1;
    
end


figure
loglog(round(logspace(0,log10(length(J)),50)),(J(round(logspace(0,log10(length(J)),50)))),'-s','linewidth',1)
xlabel('$t$','interpreter','latex')
ylabel('$\|w(t)-\theta_p\|/n$','interpreter','latex')

