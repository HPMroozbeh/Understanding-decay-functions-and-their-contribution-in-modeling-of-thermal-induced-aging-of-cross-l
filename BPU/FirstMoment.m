function y=FirstMoment(n1,n2,Rbar,q)

P=@(n)n.*((3./(2*n)).^1.5.*q.^n.*exp(-1.5*(Rbar^2./n)));
y=integral(P,n1,n2);