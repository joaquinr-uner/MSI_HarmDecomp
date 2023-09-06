function [s_imp] = impute_lsw(s,st,L,params)
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
else
    if isfield(params,'pn')
        pn = params.pn;
    else
        pn = mfilename('fullpath');
        pn = fullfile(pn(1:end-11), 'aux-functs');
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

    sp = s_imp(1:inti(1)-1);

    Left = max([1,st(qi)-3*L(qi)]);
    si = s(Left:st(qi));
    %si = sp;
    Ni = length(si);
    Li = L(qi);

    
    %saveR('temp/MissingSamples.R','si','T','sti','Li','Ni')
    name = tempname(pn_temp);
    saveR([name '.R'],'si','Li','Ni')
    %CurrentDirectory=strrep(pwd,'\','/');

   evalc(['!/usr/bin/Rscript ' pn '/RunLSW.R ' name '.R']);
   %evalc(['!"C:\Program Files\R\R-3.6.3\bin\Rscript.exe" ' pn '/RunLSW.R ' name '.R']);

    load([name '.R.mat'],'imp')
        
%    p = prctile(imp,95);
% 
%     for j=1:Li
%         if imp(j)>p
%             if j<Li 
%                 if j>1
%                     imp(j) = (imp(j-1)+imp(j+1))/2;
%                 else
%                     imp(j) = imp(j+1);
%                 end
%             else
%                 imp(j) = imp(j-1);
%             end
%         end
%     end
    s_imp(inti) = imp - mean(imp);
    intp = inti;
end

%warning('off')
%rmdir(pn_temp,'s')
%warning('on')