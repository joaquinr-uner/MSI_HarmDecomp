function [s_ext,spbk,spfw] = extendSig(s,phi,c,Np,method,opoptions)
%%
% signal extension by ARIMA model forecasting
% Inputs: 
%         s: signal to extend
%         phi: fundamental phase function of s
%         c: number of cycles for ARIMA model fitting
%         Np: number of samples to forecast
%         method: choose between forward ('fw') or backward ('bw')
%         forecasting
%         opoptions: additional options for the optimization method
% Outputs:
%         s_ext: extended signal
%         spbk: backward extended segment
%         spfw: forward extended segment

if nargin<6
    opoptions = optimoptions('fmincon');
end

s = s(:);
phi = phi(:);

N = length(s);
wrapphi = wrapTo2Pi(2*pi*phi);
dwphi = diff(wrapphi);
plocs = find(dwphi<-4.5);
pw = diff(plocs);

%[~,plocs] = findpeaks(wrapTo2Pi(2*pi*phi));
%plocs = [1; plocs; N];

msplt = strsplit(method,'-');
spfw = [];
spbk = [];
if sum(ismember(msplt,'fw'))
    %indfw = find(phi>phi(end)-c,1);
    %estInd = (indfw:N)';
    if length(pw)<c
        Seas1 = floor(median(pw));
    else
        Seas1 = floor(median(pw(end-c+1:end)));
    end
    if floor(1.1*c*Seas1)>N
        estInd = [1:N]';
    else
        estInd = (N-floor(1.1*c*Seas1):N)';
    end

    Mdlfw = regARIMA('D',0,'Seasonality',Seas1,'MALags',c,'SMALags',Seas1,'Intercept',0);
    Mdlfwest = estimate(Mdlfw,s(estInd),'X',estInd,'Display','off','Options',opoptions);
    spfw = forecast(Mdlfwest,Np,'X0',(estInd),'Y0',s(estInd),'XF',(N+1:N+Np)');
    spfw = detrend(spfw);
end

if sum(ismember(msplt,'bw'))
    %indbk = find(phi>phi(1)+c,1);
    %estIndbk = (N-indbk+1:N)';
    if length(pw)<c
        Seas2 = floor(median(pw));
    else
        Seas2 = floor(median(pw(1:c)));
    end
    sbk = flipud(s);

    estIndbk = (N-floor(1.1*c*Seas2):N)';

    Mdlbk = regARIMA('D',0,'Seasonality',Seas2,'MALags',c,'SMALags',Seas2,'Intercept',0);
    Mdlbkest = estimate(Mdlbk,sbk(estIndbk),'X',estIndbk,'Display','off','Options',opoptions);
    spbk = forecast(Mdlbkest,Np,'X0',(estIndbk),'Y0',sbk(estIndbk),'XF',(N+1:N+Np)');
    spbk = detrend(spbk);
end
    s_ext = [flipud(spbk); s; spfw];
    
end