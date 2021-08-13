function y=NIntegral(n1,n2,Rbar,q,Lambda_D)

P=@(n)((3./(2*n)).^1.5.*q.^n.*exp(-1.5*(Rbar^2./n))).*((3.*((Lambda_D.*Rbar)./(n)))./(1-(((Lambda_D.*Rbar)./(n))).^3));
y=integral(P,n1,n2);