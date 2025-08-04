function [s_imp] = impute_ddtfa(s,st,L,params)
%%
% Missing data imputation by ARIMA model forward forecasting.
% Inputs: 
%         s: signal with missing data.
%         st: initial indexes of missing data intervals.
%         L: lengths of missing data intervals.
%         params: optional parameters struct. Includes:
%                 cycl: number of cycles used to fit ARIMA model
%                 fmax: maximum STFT frequency, as proportion of sampling
%                 frequency.
%                 sigma: window parameter for STFT
% Outputs:
%         s_imp: signal with imputed values.
if nargin<4
    [A,~] = harm_decomp(s,struct('deshape',1));
    D = size(A,1);
    [~,fh] = compute_sigma(s);
    fs = length(s);
else
    if isfield(params,'D')
        D = params.D;
    else
        [A,~] = harm_decomp(s);
        D = size(A,1);
    end
    if isfield(params,'fs')
        fs = params.fs;
    else
        fs = length(s);
    end
    if isfield(params,'fh')
        fh = params.fh;
    else
        [~,fh] = compute_sigma(s,1);
        fh = fh*fs/length(s);
    end
end
s = s(:);
s_imp = s;

N = length(s);

Ni = length(st);


t = 0:1/fs:(N-1)/fs;
%t = [0:N-1]/N;
intp = 1;
ts = [];
for qi=1:Ni
    inti = st(qi):st(qi)+L(qi)-1;
    ts = [ts intp:inti(1)-1];
    intp = inti(end)+1;
end
ts = [ts inti(end)+1:N];

skp = s(ts);

IMF = zeros(N,D);
for k=1:D
    [~,IMF(:,k)]=Decompose_MP_sparse(skp,k*fh*2*pi*t',t(ts),t);

    skp = skp - IMF(ts,k);
end

s_tfa = sum(IMF,2);

for qi=1:Ni
    s_imp(st(qi):st(qi)+L(qi)-1) = s_tfa(st(qi):st(qi)+L(qi)-1);
end    

end
