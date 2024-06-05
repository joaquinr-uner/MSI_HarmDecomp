function [AD,phiD,Trend,TFdata] = harm_decomp(x,params)
%%
% Compute the harmonic decomposition of signal x.
% Inputs: 
%         x: signal under analysis
%         fs: sampling frequency of x
%         params: optional parameters struct. Includes:
%                 'sigma': window parameter sigma.
%                 'b': half-width for harmonic reconstruction
%                 'fmax': maximum frequency for the computation of the
%                 STFT, a proportion of fs
%                 'K': number of oscillatory components
%                 'r_max': maximum admissible number of harmonic
%                 'r_opt': chosen order for harmonic decomposition
%                 'Criteria': string or cell with criteria for automatic
%                 number of harmonic estimation.
%                 'crit_params': Additional parameters for trigonometric
%                 regression criteria:
%                       When 'Criteria' == 'Wang': vector of penalization
%                                                  parameters c.
%                       When 'Criteria' == 'Kavalieris': maximum order H.
%                       When 'Criteria' == 'Rl': noise variance estimate
%                       sigma.
%                'with_trend': logical variable to indicate if the signal trend
%                has to be recovered. 
if nargin<2
    sigma = compute_sigma(x,1);
    %b = round(sqrt(log(80)*sigma/pi^2)*length(x)) + 1;
    b = round(3/pi*sqrt(sigma/2)*length(x));
    fmax = 0.5;
    crit = {'Wang'};
    crit_params = struct('c',[4,6,8,10,12]);
    with_trend = 0;
    Rmax = 50;
    K = 1;
    r_opt = 0;
else
    if isfield(params,'sigma')
        sigma = params.sigma;
    else
        sigma = compute_sigma(x,1);
    end
    if isfield(params,'b')
        b = params.b;
    else
        %b = round(sqrt(log(80)*sigma/pi^2)*length(x)) + 1;
        b = round(3/pi*sqrt(sigma/2)*length(x));
    end
    if isfield(params,'fmax')
        fmax = params.fmax;
    else
        fmax = 0.5;
    end
    if isfield(params,'with_trend')
        with_trend = params.with_trend;
    else
        with_trend = 0;
    end
    if isfield(params,'deshape')
        deshape = params.deshape;
    else
        deshape = 0;
    end
    if isfield(params,'Rmax')
        Rmax = params.Rmax;
    else
        Rmax = 50;
    end
    if isfield(params,'Criteria')
        if isa(params.Criteria,'char')
            crit = {params.Criteria};
        else
           crit = params.Criteria;
        end
        crit_params = params.crit_params;
    else
        crit = {'Wang'};
        crit_params = struct('c',[4,6,8,10,12]);
    end
    if isfield(params,'K')
        K = params.K;
    else
        K = 1;
    end
    if isfield(params,'r_opt')
        r_opt = params.r_opt;
    else
        r_opt = 0;
    end
end

N = length(x);

C = zeros(K,N);

AD = [];
phiD = [];

if sigma==0
    sigma = compute_sigma(x,1);
end

[F,sF] = STFT_Gauss(x,N,sigma,fmax);
f = 0:1/N:fmax-1/N;
U = istct_fast(F,f,0.3,0.3);

W = F.*U;
c = ridge_ext(W,0.1,0.1,10,10);
cbl = max([ones(1,N);c-b]);
cbu = min([N*ones(1,N);c+b]);
W(cbl:cbu,:) = 0;
C(1,:) = ridge_correct(c,F,b,1);

for k=2:K
    ck = ridge_ext(W(10:end,:),0.1,0.1,10,10);
    ck = ck + 10;

    C(k,:) = ridge_correct(ck,F,b,1);
    if abs(median(C(k,:))-median(C(1,:)))<2*b % Check if detected ridge is a residual of c_1
        K = K - 1;
        C(k,:) = [];
    end
end

cmin = inf;
km = 1;
for k=1:K
    mink = min(C(k,:));
    if mink<cmin
        cmin = mink;
        km = k;
    end
end

if with_trend
    ckb = max([ones(1,N);C(km,:)-2*b]);
    Trend = sqrt(2)*2/max(sF)*real(sum(F(1:ckb,:),1));
    Trend = Trend(:);
    x = x - Trend;

    [F,sF] = STFT_Gauss(x,N,sigma,fmax);
else
    Trend = zeros(N,1);
end


A = zeros(K,N);
phi = zeros(K,N);
r_max = zeros(K,1);
for k=1:K
    r_max(k) = floor(0.5*N/max(C(k,:)));
    [A(k,:),phi(k,:)] = extract_harmonics(F,sF,C(k,:),b,b,1);
end
r_max(r_max>Rmax) = Rmax;

if (max(C(:))<=0.5*N && sum(~isnan(C(:))))
    if r_opt == 0
        r_opt = order_optK(x,r_max,A,phi,crit,crit_params);
    end

    for k=1:k
        ck = C(k,:);
        [ADk,phiDk] = extract_harmonics(F,sF,ck,b,b,r_opt(k));
        AD = [AD; ADk];
        phiD = [phiD; phiDk];
    end
end

TFdata.F = F;
TFdata.sigma = sigma;
TFdata.fmax = fmax;
TFdata.c = C;
end