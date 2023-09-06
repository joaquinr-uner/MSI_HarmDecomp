function [s_int, trend_int, A_int, phi_int, CONF] = harm_int_ecg(A,phi,st,L,intp,s_imp,trend)
%%
% Missing data imputation refinement by harmonic interpolation.
% Inputs: 
%         A: harmonic amplitudes of the imputed signal.
%         phi: harmonic phases of the imputed signal.
%         st: initial indexes of missing data intervals.
%         L: lengths of missing data intervals.
%         intp: interpolation type. spline or pchip.
%         s_imp: signal with imputed values
% Outputs:
%         s_int: signal with improved missing data imputation.
%         A_int: interpolated harmonic amplitudes.
%         phi_int: interpolated harmonic phases.

if nargin< 7
    trend = zeros(size(s_imp));
end

s_imp = s_imp(:);
s_int = s_imp;
trend_int = trend;

N = length(s_imp);
Ni = length(L);

r = size(A,1);

ed = st + L - 1;
CONFAl = cell(1,Ni);
CONFTl = cell(1,Ni);
CONFPl = cell(1,Ni);
CONFAu = cell(1,Ni);
CONFPu = cell(1,Ni);
CONFTu = cell(1,Ni);
for i=1:Ni
    s_pat = 0;
    Li = L(i);
    sti = st(i);
    edi = sti + Li - 1;

    A(:,sti:edi) = 0;
    phi(:,sti:edi) = 0;
    s_int(sti:edi) = 0;
    trend_int(sti:edi) = 0;
    if i==1
        lint = 1:st(i)-1;
    else
        lint = ed(i-1)+1:st(i)-1;
    end

    if i==Ni
        rint = ed(i)+1:N;
    else
        rint = ed(i)+1:st(i+1)-1;
    end

    xq = [lint rint];
    CONFAil = zeros(r,Li);
    CONFPil = zeros(r,Li);
    CONFTil = zeros(r,Li);
    CONFAiu = zeros(r,Li);
    CONFPiu = zeros(r,Li);
    CONFTiu = zeros(r,Li);
    for j=1:r
        Aq = A(j,xq);
        phiq = phi(j,xq);

        if strcmp(intp,'gpr')
            c = 5;
            wrapphi = wrapTo2Pi(2*pi*phi(1,lint));
            dwphi = diff(wrapphi);
            plocs = find(dwphi<-4);
            pw = diff(plocs);
            if length(pw)<c
                Seas1 = floor(median(pw));
            else
                Seas1 = floor(median(pw(end-c+1:end)));
            end
            if floor(1.1*c*Seas1)>length(lint)
                estInd = lint';
            else
                estInd = (lint(end)-floor(1.1*c*Seas1):lint(end))';
            end

            gprA = fitrgp(estInd,A(j,estInd));
            gprphi = fitrgp(estInd,phi(j,estInd));
            %gprA = fitrgp(xq',A(j,xq));
            %gprphi = fitrgp(xq',phi(j,xq));
            [Ai,~,confAi] = predict(gprA,[sti:edi]');
            [phii,~,confphii] = predict(gprphi,[sti:edi]');
            A(j,sti:edi) = Ai;
            phi(j,sti:edi) = phii;
        else
            Ai = interp1(xq,Aq,1:N,intp);
            phii = interp1(xq,phiq,1:N,intp);
            A(j,sti:edi) = Ai(sti:edi);
            phi(j,sti:edi) = phii(sti:edi);
        end
        s_pat = s_pat + 2*A(j,sti:edi)'.*cos(2*pi*phi(j,sti:edi))';
        if strcmp(intp,'gpr')
            CONFAil(j,:) = confAi(:,1);
            CONFPil(j,:) = confphii(:,1);
            CONFAiu(j,:) = confAi(:,2);
            CONFPiu(j,:) = confphii(:,2);
        end
    end
    trendq = trend(xq);
        if sum(abs(trend))== 0
            trend_int = trend;
        else
            if strcmp(intp,'gpr')
            gprT = fitrgp(estInd,trend(estInd));
            [Ti,~,confTi] = predict(gprT,[st(i):ed(i)]');
            trend_int(sti:edi) = Ti;
            CONFTil = confTi(:,1);
            CONFTiu = confTi(:,2);
            else
            trendi = interp1(xq,trendq,1:N,intp);
            trend_int(sti:edi) = trendi(sti:edi);
            end
        end
    %trend_int(sti:edi) = trend(sti:edi);
    %s_int(sti:edi) = s_pat;
    s_pat = s_pat + trend_int(sti:edi);
    llen = ceil(Li/2);
    rlen = floor(Li/2);
    Dist = inf ;
    pivot = -1 ;
    for ppp = [1: sti-llen edi+1: N-llen-rlen]
        patternL = [ppp: ppp+llen-1] ;
        patternR = [ppp+llen: ppp + llen + rlen-1] ;
        s0 = s_imp([patternL patternR]) ;

        % find the most similar pattern from the remainig signal for
        % imputation
        dd = norm(s0 - s_pat) ;
        if dd < Dist ;
            Dist = dd ; pivot = ppp ;
            FIT = s_imp([patternL patternR]) ;
        end
    end

    s_int(sti:edi) = FIT;
    CONFAl{i} = CONFAil;
    CONFPl{i} = CONFPil;
    CONFTl{i} = CONFTil;
    CONFAu{i} = CONFAiu;
    CONFPu{i} = CONFPiu;
    CONFTu{i} = CONFTiu;
end
    A_int = A;
    phi_int = phi;
    CONF = {CONFAl, CONFAu, CONFPl, CONFPu, CONFTl, CONFTu};
end