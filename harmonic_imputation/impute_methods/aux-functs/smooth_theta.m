function dtheta_s=smooth_theta(f_fit,theta,alpha)
N_t=length(f_fit);
fft_fit=fft(f_fit);
M=round((theta(end)-theta(1))/(2*pi));
xf_t=(theta(end)-theta(1));
if mod(N_t,2)==0
    k_t=[0:N_t/2-1  -N_t/2:-1]'*2*pi/xf_t;
else
    k_t=[0:(N_t+1)/2-1  -(N_t-1)/2:-1]'*2*pi/xf_t;
end
km=2*pi/xf_t*(M);
fft_a1_fit=fft_fit.*cutoff_a1(abs(abs(k_t)),km*alpha);
dtheta_s=real(ifft(fft_a1_fit));

