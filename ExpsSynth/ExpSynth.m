addpath(genpath(fullfile('..','harmonic_imputation')))
N = 4000;
sigma = 1e-4;
fmax = 0.5;
redun = 1;
b = round(3/pi*sqrt(sigma/2)*N);

ratio = [0.05:0.025:0.2];
J = 100;

ImpMethods = {'TLM','LSE','DMD','GPR','ARIMAF','ARIMAB','TBATS','EDMD','LSW','DRAGO','DDTFA'};
ImpNames = ['TLM';'LSE';'DMD';'GPR';'ARF';'ARB';'TBT';'EDD';'LSW';'DGO';'TFA'];
NM = length(ImpMethods);
params_int = struct('sigma',sigma,'b',b,...
    'fmax',fmax,'Criteria',{'Wang'},...
    'crit_params',[4,6,8,12],'redun',redun);
Ni = 3;
p_tlm = struct('k',0.5);
p_arimaf = struct('cycl',3,'fmax',fmax,'redun',redun);
p_arimab = struct('cycl',3,'fmax',fmax,'redun',redun);
p_lse = struct();
p_dmd = struct();
p_gpr = struct();
p_tbats = struct();
p_edmd = struct();
p_lsw = struct();
p_drago = struct();
p_ddtfa = struct();
for i=1:length(ratio)
    L = round(ratio(i)*N);
    mse_imp = zeros(J,NM);
    TRUE = zeros(J,N);
    S_MS = zeros(J,N);
    S_Imp = zeros(J,N);
    S_Spl = zeros(J,N);
    S_Pch = zeros(J,N);
    S_Lin = zeros(J,N);

    VL = zeros(J,Ni);
    St = zeros(J,Ni);
    mse_best = zeros(J,Ni);
    mse_spl = zeros(J,Ni);
    mse_pch = zeros(J,Ni);
    mse_lin = zeros(J,Ni);
    mae_best = zeros(J,Ni);
    mae_spl = zeros(J,Ni);
    mae_pch = zeros(J,Ni);
    mae_lin = zeros(J,Ni);
    rmse_best = zeros(J,Ni);
    rmse_spl = zeros(J,Ni);
    rmse_pch = zeros(J,Ni);
    rmse_lin = zeros(J,Ni);
    
    MSE_Best = zeros(J,1);
    MSE_Spl = zeros(J,1);
    MSE_Pch = zeros(J,1);
    
    bestM = zeros(J,1);

    params_imp = {p_tlm,p_lse,p_dmd,p_gpr,p_arimaf,p_arimab,p_tbats,p_edmd,p_lsw,p_drago,p_ddtfa};
    r_opt = zeros(J,NM);
    fprintf('Running for L = %i \n',L)
    for j=1:J
        K = 3 + randi(3);
        fs = 4000;
        T = N/fs;
        sigma = 1e-4;
        fmax = 0.5;
        b = round(3/pi*sqrt(sigma/2)*N);
        redun = 1;

        t = 0:1/fs:(N-1)/fs;
        f = 0:fs/N:fs*fmax-fs/N;
        % Deterministic Amplitude A(t) = sqrt(t+1)
        A = sqrt(t+1);
        
        %%
        % Random phase phi(t) = phi_det(t) + phi_rand(t) - phi_rand is a filtered
        % random gaussian process
        ff = cumsum(randn(1,N));
        IFr = ff./max(abs(ff));
        IF = smooth(IFr,240);
        phir = (cumsum(IF))/N;
        phi = 50*t + 5/(2*pi)*cos(2*pi*t) + 10*phir';
        e = 1:K;c = zeros(1,K);

        for k=2:K
            e(k) = e(k) + 0.01 - 0.02*rand(1,1);
        end
        x_c = cos(2*pi*phi);
        ALP = zeros(K-1,N);
        subs = {'Const','Poly','Bump','Tanh'};
        for k=2:K
            alpi = sample_haf(t,subs);
            x_c = x_c + alpi.*(cos(2*pi*e(k)*phi)+c(k)*sin(2*pi*e(k)*phi));
            ALP(k-1,:) = alpi;
        end

        true = A.*x_c;

        true = true - mean(true);

        true = true';
        SNR = Inf;
        r = 10^(-SNR/20)*std(x_c)*randn(N,1);
        Lmin = round(L/4);
        Lmax = round(L/2);

        st1 = round(N/4-0.05*N) + randi(0.1*N);
        st2 = round(N/2-0.05*N) + randi(0.1*N);
        st3 = round(3*N/4-0.05*N) + randi(0.1*N);

        N = length(true);
        s = true;
        
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
        s_imp = impute(s,sth,Lh,ImpMethods,params_imp);
        warning('on','all')
        mse_impj = zeros(1,NM);

        for l=1:size(s_imp,1)
            si = s_imp(l,:)';
            [mse_impj(l)] = compute_errors(true,si,sth,Lh);
        end
        [~,bestMj] = min(mse_impj);

        fprintf(['Running harmonic decomposition and interpolation on ' ImpMethods{bestMj} ' result ...\n'])
        si = s_imp(bestMj,:)';
        mse_best = compute_errors(true,si,sth,Lh);
        [AD,phiD] = harm_decomp(si);
        r_optj = size(AD,1);
        [si_spl,A_spl,phi_spl] = harm_int(AD,phiD,sth,Lh,'spline',si);
        mse_spl = compute_errors(true,si_spl,sth,Lh);

        [si_pch,A_pch,phi_pch] = harm_int(AD,phiD,sth,Lh,'pchip',si);
        mse_pch = compute_errors(true,si_pch,sth,Lh);

        fprintf(['harmonic decomposition and interpolation completed.... MSE_I: ' num2str(mse_best) '. MSE_S: ' num2str(mse_spl) '. MSE_P: ' num2str(mse_pch) '\n'])

        VL(j,:) = Lh;
        St(j,:) = sth;
        TRUE(j,:) = true;
        S_MS(j,:) = s;
        S_Imp(j,:) = si;
        S_Spl(j,:) = si_spl;
        S_Pch(j,:) = si_pch;

        r_opt(j) = r_optj;
        MSE_Best(j) = mse_best;
        MSE_Spl(j) = mse_spl;
        MSE_Pch(j) = mse_pch;
        
    end

    S = struct('MSE_Best',MSE_Best,'MSE_Spl',MSE_Spl,'MSE_Pch',MSE_Pch,...
                'True',TRUE,'S_MS',S_MS,'VL',VL,'St',St,'S_Imp',S_Imp,...
                'S_Spl',S_Spl,'S_Pch',S_Pch,'BestM',bestM,'ImpNames',ImpNames);
    save(['Results_MissingDataMultiInt_' num2str(L) '_phival' num2str(1) '.mat'],'S')
end
