addpath(genpath('/home/sentey/Dropbox/Github'))
addpath(genpath('/home/jruiz'))
drt = '/home/sentey/Documentos/Missing Data Imputation - Papers y Datos/TIDIS/Señales Separadas/NasalPressure';
drt_r = '/media/Datos/joaquinruiz/MissingDataReal/NasalPressure3';

severity = 'Normal';
files = dir([drt '/' severity]);

startindex = readmatrix(fullfile(drt,severity,'startindex.csv'));

files = files(3:end-1);
J = length(files);
opoptions = optimoptions(@fmincon,'Algorithm','interior-point','MaxFunctionEvaluations',100,'MaxIterations',30);
%opoptions = optimoptions('fmincon');
ratio = [0.05:0.05:0.2];
ImpMethods = {'TLM','LSE','DMD','GPR','ARIMAF','ARIMAB','TBATS'};
ImpNames = ['TLM';'LSE';'DMD';'GPR';'ARF';'ARB';'TBT'];
ErrorCrit = {'mae','mse','rmse'};
ErrorNames = ['mae';'mse';'rme'];
NE = size(ErrorCrit,2);
%ImpMethods = {'TLM','LSE','DMD','GPR'};
%ImpNames = ['TLM','LSE','DMD','GPR'];
NM = length(ImpMethods);
I = length(ratio);
%ratio = 0.2;

for j=1:J
    name = files(j).name;
    S = load(fullfile(drt,severity,name));
    fs = S.S.fs;
    T = 120;
    N = T*fs;
    sigma = 5e-7;
    fmax = 0.1;
    b = round(3/pi*sqrt(sigma/2)*N);
    redun = 1;
    rmax = 50;
    
    ini = startindex(j);
    s = S.S.signal(ini:ini+N-1);
    s = s - mean(s);

    fprintf(['Processing ' name '. Start Index: ' num2str(ini) '\n'])
    t = 0:1/fs:(N-1)/fs;
    f = 0:fs/N:fs*fmax-fs/N;

    tr = s';

    [~,fh] = compute_sigma(s);

    Th = N/fh;
    Mi = floor(0.9*Th);
    Ki = ceil(3*Th);


    %params_decomp = struct('fmax',fmax,'Criteria',{'Wang'},...
    %    'crit_params',[4,6,8,12],'r_max',rmax,'redun',redun);
    params_decomp = struct('r_max',rmax,'fmax',fmax,'with_trend',1,'deshape',0);
    Ni = 3;
    p_tlm = struct();
    M = 25;
    p_lse = struct('M',Mi,'K',Ki);
    p_dmd = struct('M',Mi,'K',Ki);
    p_gpr = struct('M',Mi,'K',Ki);
    p_arimaf = struct('cycl',3,'fmax',fmax,'redun',redun,'options',opoptions);
    p_arimab = struct('cycl',3,'fmax',fmax,'redun',redun,'options',opoptions);
    p_tbats = struct('pn','/home/sentey/Dropbox/Github/harmonic_imputation/impute_methods/aux-functs');

    params_imp = {p_tlm,p_lse,p_dmd,p_gpr,p_arimaf,p_arimab,p_tbats};

    N = length(tr);
    TRUE = zeros(I,N);
    S_MS = zeros(I,N);
    S_Imp = zeros(I,NM,N);
    S_Spl = zeros(I,NM,N);
    S_Pch = zeros(I,NM,N);
    S_Lin = zeros(I,NM,N);

    r_opt = zeros(I,NM);
    VL = zeros(I,Ni);
    St = zeros(I,Ni);
    Err_Imp = zeros(I,NM,NE);
    Err_Spl = zeros(I,NM,NE);
    Err_Pch = zeros(I,NM,NE);
    Err_Lin = zeros(I,NM,NE);

    Times_Imp = zeros(I,NM);
    Times_Decomp = zeros(I,NM);
    Times_Spl = zeros(I,NM);
    Times_Pch = zeros(I,NM);
    Times_Lin = zeros(I,NM);

    BestI = zeros(I,1);
    BestS = zeros(I,1);
    BestP = zeros(I,1);
    BestL = zeros(I,1);
    VLh = zeros(I,Ni);
    VSth = zeros(I,Ni);
    Ropt = zeros(I,1);
    parfor i=1:I
        L = round(ratio(i)*N);
        fprintf('Running for L = %i \n',L)
        Lmin = round(L/4);
        Lmax = round(L/2);

        st1 = round(N/4-0.05*N) + randi(0.05*N);
        st2 = round(N/2-0.05*N) + randi(0.05*N);
        st3 = round(3*N/4-0.05*N) + randi(0.05*N);

        Ll = floor(randfixedsum(3,1,L,Lmin,Lmax)');
        ed1 = st1 + Ll(1) - 1;
        ed2 = st2 + Ll(2) - 1;
        ed3 = st3 + Ll(3) - 1;

        s = tr;
        s(st1:ed1) = 0;
        s(st2:ed2) = 0;
        s(st3:ed3) = 0;
        fprintf('Interv 1: %i-%i. Interv2: %i-%i. Interv3: %i-%i  \n',st1,ed1,st2,ed2,st3,ed3)
        [sth,Lh] = missing_ints(s,struct('c','x','d',0.01*fs,'t',0));
        edh = sth + Lh - 1;
        fprintf('Starting imputation step...\n')
        %s_imp = impute(s,sth,Lh,ImpMethods,params_imp);
        warning('off','all')
        [s_imp,times_impj] = impute(s,sth,Lh,ImpMethods,params_imp,1);
        warning('on','all')

        Sj_Imp = zeros(NM,N);Sj_Spl = Sj_Imp;Sj_Pch = Sj_Imp;Sj_Lin = Sj_Imp;

        R_opt = zeros(NM,1);

        Errj_Imp = zeros(NM,NE);
        Errj_Spl = zeros(NM,NE);
        Errj_Pch = zeros(NM,NE);
        Errj_Lin = zeros(NM,NE);

        times_hdecomp = zeros(NM,1);times_splj = zeros(NM,1);
        times_pchj = zeros(NM,1);times_linj = zeros(NM,1);
        for k=1:NM
            fprintf(['Running harmonic decomposition and interpolation on ' ImpMethods{k} ' result ...\n'])
            sk = s_imp(k,:)';
            Errj_Imp(k,:) = compute_errors(tr,sk,sth,Lh,ErrorCrit);

            tic
            [ADk,phiDk,trendk] = harm_decomp(sk,params_decomp);
            r_optjk = size(ADk,1);
            t_hdecomp = toc;

            if isnan(ADk)
                sk_spl = zeros(1,N);
                sk_pch = zeros(1,N);
                sk_lin = zeros(1,N);

                Errj_Spl(k,:) = NaN;
                Errj_Pch(k,:) = NaN;
                Errj_Lin(k,:) = NaN;

                t_hdecomp = NaN;
                t_spl = NaN;
                t_pch = NaN;
                t_lin = NaN;
            else
                tic
                sk_spl = harm_int(ADk,phiDk,sth,Lh,'spline',sk,trendk);
                t_spl = toc;
                Errj_Spl(k,:) = compute_errors(tr,sk_spl,sth,Lh,ErrorCrit);
                tic
                sk_pch = harm_int(ADk,phiDk,sth,Lh,'pchip',sk,trendk);
                t_pch = toc;
                Errj_Pch(k,:) = compute_errors(tr,sk_pch,sth,Lh,ErrorCrit);

                tic
                sk_lin = harm_int(ADk,phiDk,sth,Lh,'lin',sk,trendk);
                t_lin = toc;
                Errj_Lin(k,:) = compute_errors(tr,sk_lin,sth,Lh,ErrorCrit);

            end
            Sj_Imp(k,:) = sk;
            Sj_Spl(k,:) = sk_spl;
            Sj_Pch(k,:) = sk_pch;
            Sj_Lin(k,:) = sk_lin;
            R_opt(k) = r_optjk;

            times_hdecomp(k) = t_hdecomp;
            times_splj(k) = t_spl;
            times_pchj(k) = t_pch;
            times_linj(k) = t_lin;
        end
        %fprintf(['harmonic decomposition and interpolation completed.... MSE_I: ' num2str(Tmse_bestj) '. MSE_S: ' num2str(Tmse_splj) '. MSE_P: ' num2str(Tmse_pchj) '\n'])

        VL(i,:) = Lh;
        St(i,:) = sth;
        TRUE(i,:) = tr;
        S_MS(i,:) = s;
        S_Imp(i,:,:) = Sj_Imp;
        S_Spl(i,:,:) = Sj_Spl;
        S_Pch(i,:,:) = Sj_Pch;
        S_Lin(i,:,:) = Sj_Lin;

        r_opt(i,:) = R_opt;

        Err_Imp(i,:,:) = Errj_Imp;
        Err_Spl(i,:,:) = Errj_Spl;
        Err_Pch(i,:,:) = Errj_Pch;
        Err_Lin(i,:,:) = Errj_Lin;

        [~,BestI(i,:)] = min(Errj_Imp(:,1));

        [~,BestS(i,:)] = min(Errj_Spl(:,1));

        [~,BestP(i,:)] = min(Errj_Pch(:,1));

        [~,BestL(i,:)] = min(Errj_Lin(:,1));

        fprintf(['Best Imputation Method: ' ImpMethods{BestI(i,:)} '. Best Imputation+Pchip_Interp Method: ' ImpMethods{BestP(i,:)} '\n'])
        Times_Decomp(i,:) = times_hdecomp;
        Times_Imp(i,:) = times_impj;
        Times_Spl(i,:) = times_splj;
        Times_Pch(i,:) = times_pchj;
        Times_Lin(i,:) = times_linj;
    end
    St = struct('Err_Imp',Err_Imp,'Err_Spl',Err_Spl,'Err_Pch',Err_Pch,'Err_Lin',...
        Err_Lin,'True',TRUE,'S_MS',S_MS,'VL',VL,'St',St,'S_Imp',S_Imp,...
        'S_Spl',S_Spl,'S_Pch',S_Pch,'S_Lin',S_Lin,'ImpNames',ImpNames, ...
        'Metrics',ErrorNames,'Times_Imp',Times_Imp,'Times_Decomp',Times_Decomp,...
        'Times_Spl',Times_Spl,'Times_Pch',Times_Pch,...
        'Times_Lin',Times_Lin,'BestI',BestI,'BestS',BestS,'BestP',BestP,...
        'BestL',BestL);
        save(fullfile(drt_r,severity,['Results_MissingDataNP_' name '.mat']),'St')
end