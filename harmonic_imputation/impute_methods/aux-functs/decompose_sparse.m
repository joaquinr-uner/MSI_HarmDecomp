function [a1_c,a1_s,a0,IMF]=decompose_sparse(f,theta,theta_f,alpha)

 [basis_theta,basis_f,w,mk,nk]=basis_theta_sparse(theta,theta_f,alpha);

 
 quiet=1;
lambda=1e-1;
pcgmaxi = 15;
[x,status]=l1_ls(basis_theta,f,lambda,1e-2,quiet,[],pcgmaxi);

 
 x0=x.*w;

a0=basis_f*x0;

xc=zeros(size(basis_theta,2),1);
xs=xc;

xc(1)=x(2*nk);
xs(1)=x(2*nk+1);
for j=1:mk
    xc(2*j)=(x(2*(nk-j))+x(2*(nk+j)));
    xc(2*j+1)=(x(2*(nk+j)+1)-x(2*(nk-j)+1));
    xs(2*j+1)=(x(2*(nk-j))-x(2*(nk+j)));
    xs(2*j)=(x(2*(nk+j)+1)+x(2*(nk-j)+1));
end

a1_c=basis_f*xc;
a1_s=basis_f*xs;

id=[(2*nk-2*mk):(2*nk+2*mk+1)];
IMF=basis_f(:,id)*x(id);