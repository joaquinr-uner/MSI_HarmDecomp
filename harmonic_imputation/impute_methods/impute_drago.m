function [s_imp] = impute_drago(s,st,L,params)
%%
% Missing data imputation by forecasting using dynamical mode 
% decomposition.
% Inputs: 
%         s: signal with missing data.
%         st: initial indexes of missing data intervals.
%         L: lengths of missing data intervals.         
%         params: optional parameters struct. Includes:
%                  K: size of the dataset for forecasting
%                  M: lengths of the segments used for forecasting
% Outputs:
%         s_imp: signal with imputed values.
s = s(:);
N = length(s);

miss_int = [];
for i=1:length(st)
    miss_int = [miss_int st(i):st(i)+L-1 ];
end

obs_int = setdiff(1:N,miss_int);
S = eye(N);
S = S(obs_int,:);

D = circulant([1 -2 1 zeros(1,N-3)]);

l1 = 0.01;
l2 = 0.02;

xi = S*s;
I = 10;
for i=1:I
    A = diag(l1./(abs(xi)+1));
    fip1 = (S'*A*S+l2*(D'*D))\S'*A*S*s;
    xip1 = abs(xi)./(abs(xi)+l1).*(S*(s-fip1));
    xi = xip1;
    s_imp = fip1;
end

end