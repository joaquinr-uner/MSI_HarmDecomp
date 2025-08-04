function [s_imp] = impute_tbats(s,st,L,params)
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
    pn = cd;
    Ti = 0;
    j = 0;
else
    if isfield(params,'pn')
        pn = params.pn;
    else
        pn = mfilename('fullpath');
        pn = [pn(1:end-13) '/aux-functs'];
    end
    if isfield(params,'T')
        T = params.T;
    else
        T = 0;
    end
end
s = s(:);

Ni = length(st);

s_imp = s;
if ~isfolder([pn '/temp'])
    mkdir([pn '/temp'])
end
pn_temp = [pn '/temp/'];
intp = 0;

for qi=1:Ni

    inti = st(qi):st(qi)+L(qi)-1;

    %sp = s_imp(1:inti(1)-1);
    sp = s(intp(end)+1:inti(1)-1);

    Np = length(sp);
    if T == 0
        [sigma,~,b] = compute_sigma(sp,1);

        ff = 0:1/length(sp):0.5-1/length(sp);
        [Ff,sFf] = STFT_Gauss(sp,Np,sigma,0.5);
        Uf = istct_fast(Ff,ff,0.3,0.03);
        Wf = Ff.*Uf;

        cf = ridge_ext(Wf,0.1,0.1,10,10);
        cf = ridge_correct(cf,Ff,b,1);

        [~,phif] = extract_harmonics(Ff,sFf,cf,b,b,1);

        wrapphi = wrapTo2Pi(2*pi*phif);dwphi = diff(wrapphi);
        plocs = find(dwphi<-5);
        pw = diff(plocs);
        if length(pw)<3
            Ti = floor(median(pw));
            ind = sum(pw);
        else
            Ti = floor(median(pw(end-2:end)));
            ind = Np-sum(pw(end-2:end)):Np;
        end
    else
        ind = Np-3*T:Np;
    end
    ind(ind<1) = [];

    si = sp(ind);
    Ni = length(si);
    sti = st(qi);
    Li = L(qi);

    %saveR('temp/MissingSamples.R','si','T','sti','Li','Ni')
    name = tempname(pn_temp);
    saveR([name '.R'],'si','Ti','sti','Li','Ni','j')
    %CurrentDirectory=strrep(pwd,'\','/');

    evalc(['!/usr/bin/Rscript ' pn '/RunTBATS.R ' name '.R']);

    load([name '.R.mat'],'imp')


    s_imp(inti) = imp - mean(imp);
    intp = inti;
end

%warning('off')
%rmdir(pn_temp,'s')
%warning('on')