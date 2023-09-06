function y=integral_fft(f)
i=sqrt(-1);
N=length(f);
ft=fft(f);
if mod(N,2)
    k=[0:N/2  -N/2:-1]'*2*pi;
else
    k=[0:N/2-1  -N/2:-1]'*2*pi;
end
fy=zeros(N,1);
fy(2:end)=-i*ft(2:end)./k(2:end);
yp=real(ifft(fy));
yl=round(ft(1)/(2*N*pi))*2*pi*[0:N]'/N;
y=[yp;yp(1)]+yl;
y=y-y(1);