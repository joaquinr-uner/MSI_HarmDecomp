function y=basis_FEM_linear(mesh_b,x)
if size(x,1)<size(x,2)
    x=x';
end
if size(mesh_b,1)<size(mesh_b,2)
    mesh_b=mesh_b';
end
n=length(x);
m=length(mesh_b);


y=zeros(n,m);
for j=1:m-1
    id=(x>=mesh_b(j)& x<=mesh_b(j+1));
    y(id,j)=-(x(id)-mesh_b(j+1))/(mesh_b(j+1)-mesh_b(j));
    y(id,j+1)=(x(id)-mesh_b(j))/(mesh_b(j+1)-mesh_b(j));
end

