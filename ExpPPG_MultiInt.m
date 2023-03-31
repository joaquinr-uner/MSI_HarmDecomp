addpath(genpath('/home/sentey/Dropbox/Github'))
addpath(genpath('/home/jruiz'))
%drt = '/home/sentey/Documentos/Missing Data Imputation - Papers y Datos/TIDIS/Señales Separadas/Plethysmogram';
%drt_r = '/media/Datos/joaquinruiz/MissingDataReal/Plethysmogram';
drt = 'Plethysmogram';
drt_r = 'Results/Plethysmogram';
severity = 'Normal';
files = dir([drt '/' severity]);

%startindex = readmatrix(fullfile(drt,severity,'startindex.csv'));
load(fullfile(drt,severity,'startindex.mat'))
files = files(3:end-1);
J = length(files);

ratio = [0.05:0.025:0.2];
ImpMethods = {'TLM','LSE','DMD','GPR','ARIMAF','ARIMAB'};
ImpNames = ['TLM';'LSE';'DMD';'GPR';'ARF';'ARB'];
NM = length(ImpMethods);
%ImpMethods = {'TLM','LSE','DMD'};
%ImpNames = ['TLM','LSE','DMD'];
I = length(ratio);
%ratio = 0.2;
for j=1:J
    name = files(j).name;
    load(fullfile(drt,severity,name))
    fs = S.fs;
    T = 120;
    N = T*fs;
    sigma = 1e-6;
    fmax = 0.1;
    b = round(3/pi*sqrt(sigma/2)*N);
    redun = 1;
    rmax = 50;

    params_int = struct('sigma',sigma,'b',b,...
        'fmax',fmax,'Criteria',{'Wang'},...
        'crit_params',[4,6,8,12],'rmax',rmax,'redun',redun);
    Ni = 3;
    p_tlm = struct();
    M = 25;
    p_lse = struct();
    p_dmd = struct();
    p_gpr = struct();
    p_arimaf = struct('cycl',3,'fmax',fmax,'redun',redun);
    p_arimab = struct('cycl',3,'fmax',fmax,'redun',redun);

    params_imp = {p_tlm,p_lse,p_dmd,p_gpr,p_arimaf,p_arimab};
    ini = startindex(j);
    s = S.signal(ini:ini+N-1);
    s = s - mean(s);

    fprintf(['Processing ' name '. Start Index: ' num2str(ini) '\n'])
    t = 0:1/fs:(N-1)/fs;
    f = 0:fs/N:fs*fmax-fs/N;

    true = s';

    N = length(true);
    MSE_Imp = zeros(I,NM);
    MSE_Best = zeros(I,Ni);
    MSE_Spl = zeros(I,Ni);
    MSE_Pch = zeros(I,Ni);
    MSE_Lin = zeros(I,Ni);
    MAE_Best = zeros(I,Ni);
    MAE_Spl = zeros(I,Ni);
    MAE_Pch = zeros(I,Ni);
    MAE_Lin = zeros(I,Ni);
    RMSE_Best = zeros(I,Ni);
    RMSE_Spl = zeros(I,Ni);
    RMSE_Pch = zeros(I,Ni);
    RMSE_Lin = zeros(I,Ni);
    TMSE_Best = zeros(I,1);
    TMSE_Spl = zeros(I,1);
    TMSE_Pch = zeros(I,1);
    TMSE_Lin = zeros(I,1);
    TMAE_Best = zeros(I,1);
    TMAE_Spl = zeros(I,1);
    TMAE_Pch = zeros(I,1);
    TMAE_Lin = zeros(I,1);
    TRMSE_Best = zeros(I,1);
    TRMSE_Spl = zeros(I,1);
    TRMSE_Pch = zeros(I,1);
    TRMSE_Lin = zeros(I,1);

    TRUE = zeros(I,N);
    Sigs = zeros(I,N);
    SI = zeros(I,N);
    S_Spl = zeros(I,N);
    S_Pch = zeros(I,N);
    S_Lin = zeros(I,N);
    BestM = zeros(I,1);
    VLh = zeros(I,Ni);
    VSth = zeros(I,Ni);
    Ropt = zeros(I,1);
    parfor i=1:I
        L = round(ratio(i)*N);
        fprintf('Running for L = %i \n',L)
        Lmin = round(L/4);
        Lmax = round(L/2);

        st1 = round(N/4-0.05*N) + randi(0.1*N);
        st2 = round(N/2-0.05*N) + randi(0.1*N);
        st3 = round(3*N/4-0.05*N) + randi(0.1*N);

        Ll = floor(randfixedsum(3,1,L,Lmin,Lmax)');
        ed1 = st1 + Ll(1) - 1;
        ed2 = st2 + Ll(2) - 1;
        ed3 = st3 + Ll(3) - 1;

        s = true;
        s(st1:ed1) = 0;
        s(st2:ed2) = 0;
        s(st3:ed3) = 0;

        fprintf('Interv 1: %i-%i. Interv2: %i-%i. Interv3: %i-%i  \n',st1,ed1,st2,ed2,st3,ed3)
        [sth,Lh] = missing_ints(s,0.01*fs,0);
        edh = sth + Lh - 1;
        fprintf('Starting imputation step...\n')
        s_imp = impute(s,sth,Lh,ImpMethods,params_imp);

        mse_imp = zeros(1,NM);

        for l=1:size(s_imp,1)
            si = s_imp(l,:)';
            [mse_imp(l)] = compute_errors(true,si,sth,Lh);
        end
        [~,bestM] = min(mse_imp);

        fprintf(['Running harmonic decomposition and interpolation on ' ImpMethods{bestM} ' result ...\n'])
        si = s_imp(bestM,:)';
        [Tmse_best, mse_best, Tmae_best, mae_best, Trmse_best, rmse_best] = compute_errors(true,si,sth,Lh);
        [AD,phiD] = harm_decomp(si,params_int);
        r_opt = size(AD,1);
        [si_spl,A_spl,phi_spl] = harm_int(AD,phiD,sth,Lh,'spline',si);
        [Tmse_spl,mse_spl, Tmae_spl,mae_spl, Trmse_spl,rmse_spl] = compute_errors(true,si_spl,sth,Lh);
        
        [si_pch,A_pch,phi_pch] = harm_int(AD,phiD,sth,Lh,'pchip',si);
        [Tmse_pch,mse_pch, Tmae_pch,mae_pch, Trmse_pch,rmse_pch] = compute_errors(true,si_pch,sth,Lh);
        
        [si_lin,A_lin,phi_lin] = harm_int(AD,phiD,sth,Lh,'lin',si);
        [Tmse_lin,mse_lin, Tmae_lin,mae_lin, Trmse_lin,rmse_lin] = compute_errors(true,si_lin,sth,Lh);


        fprintf('harmonic decomposition and interpolation completed...\n')

        MSE_Imp(i,:) = mse_imp;
        MSE_Best(i,:) = mse_best;
        MSE_Spl(i,:) = mse_spl;
        MSE_Pch(i,:) = mse_pch;
        MSE_Lin(i,:) = mse_lin;
        MAE_Best(i,:) = mae_best;
        MAE_Spl(i,:) = mae_spl;
        MAE_Pch(i,:) = mae_pch;
        MAE_Lin(i,:) = mae_lin;
        RMSE_Best(i,:) = rmse_best;
        RMSE_Spl(i,:) = rmse_spl;
        RMSE_Pch(i,:) = rmse_pch;
        RMSE_Lin(i,:) = rmse_lin;
        TMSE_Best(i) = Tmse_best;
        TMSE_Spl(i) = Tmse_spl;
        TMSE_Pch(i) = Tmse_pch;
        TMSE_Lin(i) = Tmse_lin;
        TMAE_Best(i) = Tmae_best;
        TMAE_Spl(i) = Tmae_spl;
        TMAE_Pch(i) = Tmae_pch;
        TMAE_Lin(i) = Tmae_lin;
        TRMSE_Best(i) = Trmse_best;
        TRMSE_Spl(i) = Trmse_spl;
        TRMSE_Pch(i) = Trmse_pch;
        TRMSE_Lin(i) = Trmse_lin;

        TRUE(i,:) = true;
        Sigs(i,:) = s;
        SI(i,:) = si;
        S_Spl(i,:) = si_spl;
        S_Pch(i,:) = si_pch;
        S_Lin(i,:) = si_lin;
        BestM(i) = bestM;
        VLh(i,:) = Lh;
        VSth(i,:) = sth;
        Ropt(i) = r_opt;
    end

        St = struct('MSE_Imp',MSE_Imp,'MSE_Best',MSE_Best,'MSE_Spl',MSE_Spl,...
            'MSE_Pch',MSE_Pch,'MSE_Lin',MSE_Lin,...
            'MAE_Best',MAE_Best,'MAE_Spl',MAE_Spl,'MAE_Pch',MAE_Pch, ...
            'MAE_Lin',MAE_Lin,'RMSE_Best',RMSE_Best,'RMSE_Spl',RMSE_Spl, ...
            'RMSE_Pch',RMSE_Pch,'RMSE_Lin',RMSE_Lin, ...
            'TMSE_Best',TMSE_Best,'TMSE_Spl',TMSE_Spl,'TMSE_Pch',TMSE_Pch, ...
            'TMSE_Lin',TMSE_Lin,'TMAE_Best',TMAE_Best,'TMAE_Spl',TMAE_Spl, ...
            'TMAE_Pch',TMAE_Pch,'TMAE_Lin',TMAE_Lin,'TRMSE_Best',TRMSE_Best, ...
            'TRMSE_Spl',TRMSE_Spl,'TRMSE_Pch',TRMSE_Pch,'TRMSE_Lin',TRMSE_Lin,...
            'True',TRUE,'S_MS',Sigs,'VL',VLh,'St',VSth,'S_Imp',SI,...
            'S_Spl',S_Spl,'S_Pch',S_Pch,'S_Lin',S_Lin,'BestM',BestM,...
            'Ropt',Ropt,'ImpNames',ImpNames);
        save(fullfile(drt_r,severity,['Results_MissingDataPPG_' name '.mat']),'St')
end