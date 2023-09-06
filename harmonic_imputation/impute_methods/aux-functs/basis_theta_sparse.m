function [y,yf,w,mk,n]=basis_theta_sparse(theta,theta_f,alpha)
p2=2*pi;

n=1*round((theta_f(end)-theta_f(1))/p2);
m=min(4*n,length(theta_f));
mk=floor(alpha*n);
%m=n;

y=zeros(length(theta),2*m+1);
yf=zeros(length(theta_f),2*m+1);

y(:,1)=1;
yf(:,1)=1;

w=zeros(2*m+1,1);
w(1)=1;

for k=1:m
    y(:,2*k)=cos(k*theta/(n));
    y(:,2*k+1)=sin(k*theta/(n));
    yf(:,2*k)=cos(k*theta_f/(n));
    yf(:,2*k+1)=sin(k*theta_f/(n));
end

if mk>0
    for j=1:mk
        w(2*j)=1;
        w(2*j+1)=1;
    end
    
end