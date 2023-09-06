function [s_ext,spbk,spfw] = extend_fw(s,phi,c,Np)
s = s(:);
phi = phi(:);

N = length(s);
[~,plocs] = findpeaks(wrapTo2Pi(2*pi*phi));
plocs = [1; plocs; N];

pw = diff(plocs);
indfw = find(phi>phi(end)-c,1);
Seas1 = round(median(pw(end-c+1:end)));
estInd = (indfw:N)';
Mdlfw = regARIMA('D',0,'Seasonality',Seas1,'MALags',c,'SMALags',Seas1,'Intercept',0);
Mdlfwest = estimate(Mdlfw,s(estInd),'X',estInd,'Display','off');
spfw = forecast(Mdlfwest,Np,'X0',(estInd),'Y0',s(estInd),'XF',(N+1:N+Np)');
spfw = detrend(spfw);

s_ext = [s; spfw];
end