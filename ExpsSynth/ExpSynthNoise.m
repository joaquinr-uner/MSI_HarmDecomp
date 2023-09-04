#addpath(genpath('/home/sentey/Dropbox/Github'))
# Add Path to harmonic_imputation toolbox
# addpath(genpath(...))

            N = 4000;
            sigma = 1e-4;
            fs = 4000;
            fmax = 0.5;
            T = N/fs;
            redun = 1;
            b = round(3/pi*sqrt(sigma/2)*N);
ratio = [0.05:0.05:0.2];
J = 100;
ImpMethods = {'TLM','LSE','DMD','GPR','ARIMAF','ARIMAB','TBATS','DDTFA','EDMD','LSW'};
ImpNames = ['TLM';'LSE';'DMD';'GPR';'ARF';'ARB';'TBT';'TFA';'EDD';'LSW'];
ErrorCrit = {'mae','mse','rmse'};
ErrorNames = ['mae';'mse';'rme'];
NE = size(ErrorCrit,2);
vSNR = [10,20, Inf];

NM = length(ImpMethods);

params_decomp = struct(...
    'fmax',fmax,'Criteria',{'Wang'},...
    'crit_params',[4,6,8,12],'redun',redun);

Ni = 3;

p_tlm = struct('k',0.5);
p_arimaf = struct('cycl',3,'fmax',0.5,'redun',1);
p_arimab = struct('cycl',3,'fmax',0.5,'redun',1);
p_edmd = struct();
p_tbats = struct('pn','/home/sentey/Dropbox/Github/harmonic_imputation/impute_methods/aux-functs');
p_lsw= struct('pn','/home/sentey/Dropbox/Github/harmonic_imputation/impute_methods/aux-functs');
p_lse = struct();
p_dmd = struct();
p_gpr = struct();
p_ddtfa = struct('fs',fs);

for l=1:length(vSNR)
    for i=1:length(ratio)
        L = round(ratio(i)*N);
        TRUE = zeros(J,N);
        TRUEN = zeros(J,N);
        S_MS = zeros(J,N);
        S_Imp = zeros(J,NM,N);
        S_Spl = zeros(J,NM,N);
        S_Pch = zeros(J,NM,N);
        S_Lin = zeros(J,NM,N);

        VL = zeros(J,Ni);
        St = zeros(J,Ni);
        Err_Imp = zeros(J,NM,NE);
        Err_Spl = zeros(J,NM,NE);
        Err_Pch = zeros(J,NM,NE);
        Err_Lin = zeros(J,NM,NE);

        Times_Imp = zeros(J,NM);
        Times_Decomp = zeros(J,NM);
        Times_Spl = zeros(J,NM);
        Times_Pch = zeros(J,NM);
        Times_Lin = zeros(J,NM);

        BestI = zeros(J,1);
        BestS = zeros(J,1);
        BestP = zeros(J,1);
        BestL = zeros(J,1);

        %M = 25;
        %p_lse = struct('M',1.5*L,'K',round(3.75*L));
        %p_dmd = struct('M',1.5*L,'K',round(3.75*L));
        %p_gpr = struct('M',1.5*L,'K',round(3.75*L),'fmode','sd','pmode','bcd');

        params_imp = {p_tlm,p_lse,p_dmd,p_gpr,p_arimaf,p_arimab,p_tbats,p_edmd,p_lsw};
        %params_imp = {p_arimaf,p_arimab};

        r_opt = zeros(J,NM);
        fprintf('Running for L = %i \n',L)
        for j=1:J

            K = 3 + randi(3);
            t = 0:1/fs:(N-1)/fs;
            f = 0:fs/N:fs*fmax-fs/N;

            A = sqrt(t+1);

            %%
            ff = cumsum(randn(1,N));
            IFr = ff./max(abs(ff));
            IF = smooth(IFr,240);
            phir = (cumsum(IF))/N;

            phi = 50*t + 5/(2*pi)*cos(2*pi*t) + 10*phir';

            e = 1:K;c = zeros(1,K);

            for k=2:K
                e(k) = e(k) + 0.01 - 0.02*rand(1,1);%c = [0,0.25,0.25,0.25,0.25,0.25]*pi;
            end
            x_c = cos(2*pi*phi);
            ALP = zeros(K-1,N);
            subs = {'Const','Poly','Bump','Tanh'};
            %subs = {'Const'};
            for k=2:K
                alpi = sample_haf(t,subs);
                x_c = x_c + alpi.*(cos(2*pi*e(k)*phi)+c(k)*sin(2*pi*e(k)*phi));
                ALP(k-1,:) = alpi;
            end

            tr = A.*x_c;
            r = 10^(-vSNR(l)/20)*std(tr)*randn(size(tr));

            trc = tr;
            tr = tr + r;

            tr = tr - mean(tr);

            tr = tr';
            Lmin = round(L/4);
            Lmax = round(L/2);

            st1 = round(N/4-0.05*N) + randi(0.1*N);
            st2 = round(N/2-0.05*N) + randi(0.1*N);
            st3 = round(3*N/4-0.05*N) + randi(0.1*N);

            N = length(tr);
            s = tr;
            
            Ll = floor(randfixedsum(3,1,L,Lmin,Lmax)');
            ed1 = st1 + Ll(1) - 1;
            ed2 = st2 + Ll(2) - 1;
            ed3 = st3 + Ll(3) - 1;

            s(st1:ed1) = 0;
            s(st2:ed2) = 0;
            s(st3:ed3) = 0;

            fprintf('Interv 1: %i-%i. Interv2: %i-%i. Interv3: %i-%i  \n',st1,ed1,st2,ed2,st3,ed3)
            [sth,Lh] = missing_ints(s,struct('c','x','d',0.01*fs,'t',0));
            edh = sth + Lh - 1;
            fprintf('Starting imputation step...\n')
            
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
                Errj_Imp(k,:) = compute_errors(trc,sk,sth,Lh,ErrorCrit);

                tic
                [ADk,phiDk] = harm_decomp(sk,params_decomp);
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
                    [sk_spl,~,ADspl,phiDspl] = harm_int(ADk,phiDk,sth,Lh,'spline',sk);
                    t_spl = toc;
                    Errj_Spl(k,:) = compute_errors(trc,sk_spl,sth,Lh,ErrorCrit);
                    tic
                    [sk_pch,~,ADpch,phiDpch] = harm_int(ADk,phiDk,sth,Lh,'pchip',sk);
                    t_pch = toc;
                    Errj_Pch(k,:) = compute_errors(trc,sk_pch,sth,Lh,ErrorCrit);

                    tic
                    [sk_lin,~,ADlin,phiDlin] = harm_int(ADk,phiDk,sth,Lh,'lin',sk);
                    t_lin = toc;
                    Errj_Lin(k,:) = compute_errors(trc,sk_lin,sth,Lh,ErrorCrit);

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

            VL(j,:) = Lh;
            St(j,:) = sth;
            TRUE(j,:) = trc;
            TRUEN(j,:) = tr;
            S_MS(j,:) = s;
            S_Imp(j,:,:) = Sj_Imp;
            S_Spl(j,:,:) = Sj_Spl;
            S_Pch(j,:,:) = Sj_Pch;
            S_Lin(j,:,:) = Sj_Lin;

            r_opt(j,:) = R_opt;

            Err_Imp(j,:,:) = Errj_Imp;
            Err_Spl(j,:,:) = Errj_Spl;
            Err_Pch(j,:,:) = Errj_Pch;
            Err_Lin(j,:,:) = Errj_Lin;

            [~,BestI(j,:)] = min(Errj_Imp(:,1));

            [~,BestS(j,:)] = min(Errj_Spl(:,1));

            [~,BestP(j,:)] = min(Errj_Pch(:,1));

            [~,BestL(j,:)] = min(Errj_Lin(:,1));

            fprintf(['Best Imputation Method: ' ImpMethods{BestI(j,:)} '. Best Imputation+Pchip_Interp Method: ' ImpMethods{BestP(j,:)} '\n'])
            Times_Imp(j,:) = times_impj;
            Times_Decomp(j,:) = times_hdecomp;
            Times_Spl(j,:) = times_splj;
            Times_Pch(j,:) = times_pchj;
            Times_Lin(j,:) = times_linj;
        end

        S = struct('Err_Imp',Err_Imp,'Err_Spl',Err_Spl,'Err_Pch',Err_Pch,'Err_Lin',...
            Err_Lin,'True',TRUE,'TrueN',TRUEN,'S_MS',S_MS,'VL',VL,'St',St,'S_Imp',S_Imp,...
            'S_Spl',S_Spl,'S_Pch',S_Pch,'S_Lin',S_Lin,'ImpNames',ImpNames, ...
            'Metrics',ErrorNames,'Times_Imp',Times_Imp,'Times_Decomp',Times_Decomp,...
            'Times_Spl',Times_Spl,'Times_Pch',Times_Pch,...
            'Times_Lin',Times_Lin,'BestI',BestI,'BestS',BestS,'BestP',BestP,...
            'BestL',BestL);
        save(['Results_MissingDataMultiInt_' num2str(L) '_' num2str(vSNR(l)) '.mat'],'S')
    end
end
