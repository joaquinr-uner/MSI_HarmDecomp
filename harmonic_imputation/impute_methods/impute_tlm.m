function [s_imp] = impute_tlm(s,st,L,params)
%%
% Missing data imputation using Takens' Lag Map algorithm.
% Inputs:
%         s: signal with missing data.
%         st: initial indexes of missing data intervals.
%         L: lengths of missing data intervals.
%         params: optional params struct. Includes:
%                   k: left and right fitting patterns length parameter.
% Outputs:
%         s_imp: signal with imputed values.
if nargin<4
    X = fft(s);
    N = length(s);
    [~,indf] = max(abs(X(50:end/2).^2));
    indf = indf + 50;
    d = round(2/indf*N);
else
    if isfield(params,'k')
        k = params.k;
        X = fft(s);
        [~,indf] = max(abs(X(50:end/2)).^2);
        indf = indf + 50;
        N = length(s);
        d = round(k/indf*N);
    else
        if isfield(params,'d')
            d = params.d;
        else
            X = fft(s);
            [~,indf] = max(abs(X(50:end/2).^2));
            indf = indf + 50;
            N = length(s);
            d = round(2/indf*N);
        end
    end
end


s = s(:);

N = length(s);
Ni = length(st);

s_imp = s;

for qi=1:Ni

    cal = st(qi): st(qi) + L(qi) - 1 ;
    if st(qi)-d < 1 
        llen = st(qi)-1 ;
    else
        llen = d ;
    end

    if st(qi) + L(qi) + d > N
        rlen = N - (st(qi)+L(qi)-1) ;
    else
        rlen = d ;
    end

    % determine the pattern to fit
    patternL = [st(qi)-llen: st(qi)-1] ;
    patternR = [st(qi) + L(qi) + 1: st(qi) + L(qi) + rlen] ;

    X = s([patternL patternR]) ;

    Dist = inf ;
    pivot = -1 ;
    for ppp = [1: st(qi)-L(qi)-llen-rlen st(qi)+L(qi)+1: N-L(qi)-rlen-llen-1]
        patternL = [ppp: ppp+llen-1] ;
        patternR = [ppp+llen + (L(qi)-1) + 1: ppp+llen + (L(qi)-1) + rlen] ;
        X0 = s([patternL patternR]) ;

        % find the most similar pattern from the remainig signal for
        % imputation
        dd = norm(X0 - X) ;
        if dd < Dist ;
            Dist = dd ; pivot = ppp ;
            FIT = s(ppp+llen: ppp+llen+L(qi)-1) ;
        end
    end

    s_imp(cal) = FIT;

end

end
