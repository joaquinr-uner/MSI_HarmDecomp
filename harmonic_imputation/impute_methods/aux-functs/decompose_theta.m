function [IMF_fit,a1_c,a1_s]=decompose_theta(f_fit,theta,alpha)
N_t=length(f_fit);

% Fourier transform
fft_fit=fft(f_fit);


% generate wave number
xf_t=(theta(end)-theta(1));
k_t=[0:N_t/2-1  -N_t/2:-1]'*2*pi/xf_t;

% determine the wavenumber where IMF concetrate around
M=round((theta(end)-theta(1))/(2*pi));
km=2*pi/xf_t*(M);

% extract the IMF from the spectral of the signal
fft_a1_fit=fft_fit.*cutoff_a1(abs(abs(k_t)-km),km*alpha);


% translate the spectral of IMF to the original point
fft_env_fit=zeros(N_t,1);
if M>1
    for j=1:M
        fft_env_fit(N_t-(M-j))=fft_a1_fit(j);
    end
    for j=M+1:min(roof(5*M),N_t)
        fft_env_fit(j-M)=fft_a1_fit(j);
    end
elseif M==1
    fft_env_fit(1)=fft_a1_fit(2);
else
    fft_env_fit(1)=fft_a1_fit(1);
end

% compute a(t) and b(t)
env=ifft(fft_env_fit);
a1_c=2*real(env);
a1_s=-2*imag(env);

IMF_fit=real(ifft(fft_a1_fit));


