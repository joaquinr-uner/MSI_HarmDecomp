function [s_int, trend_int, A_int, phi_int] = harm_int(A,phi,st,L,intp,s_imp,trend)
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
for i=1:Ni
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
    for j=1:r
        Aq = A(j,xq);
        phiq = phi(j,xq);


        Ai = interp1(xq,Aq,1:N,intp);
        phii = interp1(xq,phiq,1:N,intp);
        A(j,sti:edi) = Ai(sti:edi);
        phi(j,sti:edi) = phii(sti:edi);
        s_int(sti:edi) = s_int(sti:edi) + 2*A(j,sti:edi)'.*cos(2*pi*phi(j,sti:edi))';
    end
    trendq = trend(xq);
    if sum(abs(trend))== 0
        trend_int = trend;
    else
        trendi = interp1(xq,trendq,[lint(1):rint(end)],intp);
        trend_int(sti:edi) = trendi(length(lint)+1:length(lint)+1+Li-1);
    end
    s_int(sti:edi) = s_int(sti:edi) + trend_int(sti:edi);
end
A_int = A;
phi_int = phi;
end
