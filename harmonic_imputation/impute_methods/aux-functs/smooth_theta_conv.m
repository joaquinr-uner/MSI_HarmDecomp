function dtheta_s=smooth_theta_conv(f_fit,theta,alpha)
N_t=length(f_fit);
fe=[f_fit;f_fit(end-1:-1:2)];
fft_fit=fft(fe);
M=2*round((theta(end)-theta(1))/(2*pi));
xf_t=2*(theta(end)-theta(1));
Ne=length(fe);
if mod(Ne,2)==0
    k_t=[0:Ne/2-1  -Ne/2:-1]'*2*pi/xf_t;
else
    k_t=[0:(Ne+1)/2-1  -(Ne-1)/2:-1]'*2*pi/xf_t;
end
km=2*pi/xf_t*(M);
fft_a1_fit=fft_fit.*cutoff_a1(abs(abs(k_t)),km*alpha);

dtheta_se=real(ifft(fft_a1_fit));
dtheta_s=dtheta_se(1:N_t);

