function [IMF_fit,a1_c,a1_s,a0_fit]=decompose_theta_meyer(f_fit,theta,alpha)
N_t=length(f_fit);
fft_fit=fft(f_fit);
M=round((theta(end)-theta(1))/(2*pi));
xf_t=(theta(end)-theta(1));

k_t=[0:N_t/2-1  -N_t/2:-1]'*2*pi/xf_t;
km=2*pi/xf_t*(M);
fft_a1_fit=fft_fit.*cutoff_a1_meyer(abs(abs(k_t)-km),km*alpha);
fft_a0_fit=fft_fit.*cutoff_a1_meyer(k_t-km,km*alpha).*cutoff_a1_meyer(-k_t-km,km*alpha);
fft_env_fit=zeros(N_t,1);
if M>1
    for j=1:M
        fft_env_fit(N_t-(M-j))=fft_a1_fit(j);
    end
    for j=M+1:min(roof(2*M),N_t)
        fft_env_fit(j-M)=fft_a1_fit(j);
    end
elseif M==1
    fft_env_fit(1)=fft_a1_fit(2);
else
    fft_env_fit(1)=fft_a1_fit(1);
end

env=ifft(fft_env_fit);
a1_c=2*real(env);
%a1_c=2*real(fft_env_fit(1))*ones(N_t,1);
a1_s=-2*imag(env);
IMF_fit=real(ifft(fft_a1_fit));
a0_fit=real(ifft(fft_a0_fit))-IMF_fit;