function [s_imp] = impute_arimab(s,st,L,params)
%%
% Missing data imputation by ARIMA model backward forecasting.
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
    cycl = 3;
    fmax = 0.5;
    redun = 1;
    opoptions = optimoptions('fmincon');
    sigma = 0;
else
    if isfield(params,'cycl')
        cycl = params.cycl;
    else
        cycl = 3;
    end
    if isfield(params,'fmax')
        fmax = params.fmax;
    else
        fmax = 0.5;
    end
    if isfield(params,'redun')
        redun = params.redun;
    else
        redun = 1;
    end
    if isfield(params,'options')
        opoptions = params.options;
    else
        opoptions = optimoptions('fmincon');
    end
    if isfield(params,'sigma')
        sigma = params.sigma;
    else
        sigma = 0;
    end
end
s = s(:);

N = length(s);
Ni = length(st);

s = flipud(s);

s_imp = s;

stf = fliplr(N-st-L+2);
Lf = fliplr(L);
intp = 0;
for qi=1:Ni

    inti = stf(qi):stf(qi)+Lf(qi)-1;

    %sp = s_imp(1:inti(1)-1);
    sp = s(intp(end)+1:inti(1)-1);
    
    if sigma==0
        sigma = compute_sigma(sp,1);
    end
    
    b = round(3/pi*sqrt(sigma/2)*length(sp));
    b = b*redun;

    ff = 0:1/length(sp):0.5-1/length(sp);
    [Ff,sFf] = STFT_Gauss(sp,length(sp)*redun,sigma,0.5);
    Uf = istct_fast(Ff,ff,0.3,0.03);
    Wf = Ff.*Uf;
    
    cf = ridge_ext(Wf,0.1,0.1,10,10);
    cf = ridge_correct(cf,Ff,b,1);

    [~,phif] = extract_harmonics(Ff,sFf,cf,b,b,1);

    %phif = phif/redun;
    Np = Lf(qi);

    ext_arima = extendSig(sp,phif,cycl,Np,'fw',opoptions);

    s_imp(inti) = ext_arima(end-Np+1:end);
    intp = inti;
end

s_imp = flipud(s_imp);
