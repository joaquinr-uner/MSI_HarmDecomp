addpath(genpath('/home/sentey/Dropbox/Github'))
addpath(genpath('/home/jruiz'))

N = 8000;
sigma = 1e-4;
fmax = 0.5;
redun = 1;
b = round(3/pi*sqrt(sigma/2)*N);

ratio = [0.05:0.025:0.2];
ratio = 0.2;
J = 100;

%ImpMethods = {'TLM','LSE','DMD','GPR','ARIMAF','ARIMAB'};
%ImpNames = ['TLM';'LSE';'DMD';'GPR';'ARF';'ARB'];
ImpMethods = {'ARIMAF','ARIMAB'};
ImpNames = ['ARF','ARB'];
NM = length(ImpMethods);
params_int = struct('sigma',sigma,'b',b,...
    'fmax',fmax,'Criteria',{'Wang'},...
    'crit_params',[4,6,8,12],'redun',redun);
Ni = 3; 
p_tlm = struct('k',0.5);
p_arimaf = struct('cycl',3,'fmax',fmax,'redun',redun);
p_arimab = struct('cycl',3,'fmax',fmax,'redun',redun);
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
    
    Tmse_best = zeros(J,1);
    Tmse_spl = zeros(J,1);
    Tmse_pch = zeros(J,1);
    Tmse_lin = zeros(J,1);
    Tmae_best = zeros(J,1);
    Tmae_spl = zeros(J,1);
    Tmae_pch = zeros(J,1);
    Tmae_lin = zeros(J,1);
    Trmse_best = zeros(J,1);
    Trmse_spl = zeros(J,1);
    Trmse_pch = zeros(J,1);
    Trmse_lin = zeros(J,1);
    
    bestM = zeros(J,1);
    M = 25;
    p_lse = struct();
    p_dmd = struct();
    p_gpr = struct();

    %p_lse = struct('M',1.5*L,'K',round(3.75*L));
    %p_dmd = struct('M',1.5*L,'K',round(3.75*L));
    %p_gpr = struct('M',1.5*L,'K',round(3.75*L),'fmode','sd','pmode','bcd');

    params_imp = {p_tlm,p_lse,p_dmd,p_gpr,p_arimaf,p_arimab};
    r_opt = zeros(J,NM);
    fprintf('Running for L = %i \n',L)
    for j=1:J
        N = 8000;
        K = 3 + randi(3);
        fs = 4000;
        T = N/fs;
        sigma = 1e-4;
        fmax = 0.5;
        b = round(3/pi*sqrt(sigma/2)*N);
        redun = 1;

        t = 0:1/fs:(N-1)/fs;
        f = 0:fs/N:fs*fmax-fs/N;
        % Random ampitude A(t) = A_det(t) + A_rand(t) - A_rand is a filtered random
        % gaussian process
        sig_A = 1e-1;
        gr = conv(randn(1,N),exp(-(t/sig_A).^2));
        Ar = cumsum(gr(1:N))/N;
        A = sqrt(t+1);
        %A = max(abs(Ar))*sqrt(t+1) + Ar;
        %A = 2*max(abs(Ar))*ones(1,N) + Ar;
        
        %%
        % Random phase phi(t) = phi_det(t) + phi_rand(t) - phi_rand is a filtered
        % random gaussian process
        sig_p = 1e-1;
        gr = conv(randn(1,N),exp(-(t/sig_p).^2));
        phir = cumsum(gr(1:N))/N;
        %phi = 50*t + 2.5/(2*pi)*cos(2*pi*t) + phir;
        %phi = 50*t + 5/(2*pi)*cos(2*pi*t);
        phi = 50*t + 5*t.^2;
        %phi = 50*t + phir;

        e = 1:K;c = zeros(1,K);

        for k=2:K
            e(k) = e(k) + 0.01 - 0.02*rand(1,1);%c = [0,0.25,0.25,0.25,0.25,0.25]*pi;
        end
        x_c = cos(2*pi*phi);
        ALP = zeros(K-1,N);
        %subs = {'Const','Poly','Bump','Tanh'};
        subs = {'Const'};
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
        [sth,Lh] = missing_ints(s,0.01*fs,0);
        edh = sth + Lh - 1;
        fprintf('Starting imputation step...\n')
        %s_imp = impute(s,sth,Lh,ImpMethods,params_imp);
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
        [Tmse_bestj, mse_bestj, Tmae_bestj, mae_bestj, Trmse_bestj, rmse_bestj] = compute_errors(true,si,sth,Lh);
        [AD,phiD] = harm_decomp(si);
        r_optj = size(AD,1);
        [si_spl,A_spl,phi_spl] = harm_int(AD,phiD,sth,Lh,'spline',si);
        [Tmse_splj,mse_splj, Tmae_splj,mae_splj, Trmse_splj,rmse_splj] = compute_errors(true,si_spl,sth,Lh);

        [si_pch,A_pch,phi_pch] = harm_int(AD,phiD,sth,Lh,'pchip',si);
        [Tmse_pchj,mse_pchj, Tmae_pchj,mae_pchj, Trmse_pchj,rmse_pchj] = compute_errors(true,si_pch,sth,Lh);

        [si_lin,A_lin,phi_lin] = harm_int(AD,phiD,sth,Lh,'lin',si);
        [Tmse_linj,mse_linj, Tmae_linj,mae_linj, Trmse_linj,rmse_linj] = compute_errors(true,si_lin,sth,Lh);

        %fprintf(['harmonic decomposition and interpolation completed.... MSE_I: ' num2str(Tmse_bestj) '. MSE_S: ' num2str(Tmse_splj) '. MSE_P: ' num2str(Tmse_pchj) '\n'])

        VL(j,:) = Lh;
        St(j,:) = sth;
        TRUE(j,:) = true;
        S_MS(j,:) = s;
        S_Imp(j,:) = si;
        S_Spl(j,:) = si_spl;
        S_Pch(j,:) = si_pch;
        S_Lin(j,:) = si_lin;

        r_opt(j) = r_optj;
        mse_imp(j,:) = mse_impj;
        mse_best(j,:) = mse_bestj;
        mse_spl(j,:) = mse_splj;
        mse_pch(j,:) = mse_pchj;
        mse_lin(j,:) = mse_linj;
        mae_best(j,:) = mae_bestj;
        mae_spl(j,:) = mae_splj;
        mae_pch(j,:) = mae_pchj;
        mae_lin(j,:) = mae_linj;
        rmse_best(j,:) = rmse_bestj;
        rmse_spl(j,:) = rmse_splj;
        rmse_pch(j,:) = rmse_pchj;
        rmse_lin(j,:) = rmse_linj;
        Tmse_best(j) = Tmse_bestj;
        Tmse_spl(j) = Tmse_splj;
        Tmse_pch(j) = Tmse_pchj;
        Tmse_lin(j) = Tmse_linj;
        Tmae_best(j) = Tmae_bestj;
        Tmae_spl(j) = Tmae_splj;
        Tmae_pch(j) = Tmae_pchj;
        Tmae_lin(j) = Tmae_linj;
        Trmse_best(j) = Trmse_bestj;
        Trmse_spl(j) = Trmse_splj;
        Trmse_pch(j) = Trmse_pchj;
        Trmse_lin(j) = Trmse_linj;
        bestM(j) = bestMj;
    end

    S = struct('MSE_Imp',mse_imp,'MSE_Best',mse_best,'MSE_Spl',mse_spl,...
                'MSE_Pch',mse_pch,'MSE_Lin',mse_lin,...
                'MAE_Best',mae_best,'MAE_Spl',mae_spl,'MAE_Pch',mae_pch, ...
                'MAE_Lin',mae_lin,'RMSE_Best',rmse_best,'RMSE_Spl',rmse_spl, ...
                'RMSE_Pch',rmse_pch,'RMSE_Lin',rmse_lin, ...
                'TMSE_Best',Tmse_best,'TMSE_Spl',Tmse_spl,'TMSE_Pch',Tmse_pch, ...
                'TMSE_Lin',Tmse_lin,'TMAE_Best',Tmae_best,'TMAE_Spl',Tmae_spl, ...
                'TMAE_Pch',Tmae_pch,'TMAE_Lin',Tmae_lin,'TRMSE_Best',Trmse_best, ...
                'TRMSE_Spl',Trmse_spl,'TRMSE_Pch',Trmse_pch,'TRMSE_Lin',Trmse_lin,...
                'True',TRUE,'S_MS',S_MS,'VL',VL,'St',St,'S_Imp',S_Imp,...
                'S_Spl',S_Spl,'S_Pch',S_Pch,'S_Lin',S_Lin,'BestM',bestM,'ImpNames',ImpNames);
    save(['Results_MissingDataMultiInt_' num2str(L) '_phival' num2str(1) '.mat'],'S')
end
