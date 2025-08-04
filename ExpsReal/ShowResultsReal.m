addpath(genpath('/home/sentey/Dropbox/Github'))
addpath(genpath('C:\Users\Intel\Dropbox\Github'))
addpath(genpath('/home/jruiz'))
severity = 'Normal';

SName = {'Plethysmogram/Normal','CHARIS','Accelerometry','Airflow/Normal','NasalPressure/Normal','Thorax/Normal'};
Labels = {'PPG','ABP','ACC','AF','NP','THO'};
K = length(SName);

ratio = 0.05:0.05:0.2;
Nl = length(ratio);
ratio_s = {'5%','10%','15%','20%'};
alpha_b = 0.05/3;

PVAL_BIP = zeros(K,Nl);
PVAL_BIS = zeros(K,Nl);
PVAL_BSP = zeros(K,Nl);
HVAL_BIP = zeros(K,Nl);
HVAL_BIS = zeros(K,Nl);
HVAL_BSP = zeros(K,Nl);
MED_IB = zeros(K,Nl);
MED_SB = zeros(K,Nl);
MED_PB = zeros(K,Nl);
YLim = [0 0.32; 0 0.37; 0 0.28; 0 0.75; 0 0.37; 0 0.40];
LabX = 0.75*ones(1,6);
LabY = [0.29, 0.335, 0.25, 0.69, 0.335, 0.36];
LabX2 = 8.5*ones(1,6);
LabY2 = [55 64 37 27.5 27.5 37];

PosB = [0.11 0.64 0.25 0.26;
        0.40 0.64 0.25 0.26;
        0.69 0.64 0.25 0.26;
        0.11 0.31 0.25 0.26;
        0.40 0.31 0.25 0.26;
        0.69 0.31 0.25 0.26;];
for k=1:K
    sname = SName{k};
    label = Labels{k};
    %drt_r = ['E:\MissingDataReal\' signal];
    drt_r = ['/media/Datos/joaquinruiz/MissingDataReal/' sname];
    files_r = dir(fullfile(drt_r));
    files_r = files_r(3:end);

    J = length(files_r);

    load(fullfile(drt_r,files_r(1).name))
    ratio = 0.05:0.05:0.2;
    Nl = length(ratio);
    ErrorCrit = upper(cellstr(St.Metrics));
    Ei = 1;
    ImpNames = cellstr(St.ImpNames);
    NM = size(ImpNames,1);
    err_imp = zeros(J,Nl,NM);
    err_spl = zeros(J,Nl,NM);
    err_pch = zeros(J,Nl,NM);

    ratio_s = {'0.5','0.1','0.15','0.2'};
    alpha_b = 0.05/3;
    pval_ip = zeros(Nl,1);pval_is = zeros(Nl,1);pval_sp = zeros(Nl,1);
    pval_il = zeros(Nl,1);pval_sl = zeros(Nl,1);pval_pl = zeros(Nl,1);
    hval_ip = zeros(Nl,1);hval_is = zeros(Nl,1);hval_sp = zeros(Nl,1);

    times_imp = zeros(Nl,NM);times_decomp = times_imp; times_pch = times_imp;
    times_spl = times_imp;

    best_err_imp = zeros(J,1); best_err_spl = best_err_imp; best_err_pch = best_err_imp;

    med_ib = zeros(Nl,1); med_sb = med_ib; med_pb = med_ib;
    pval_bip = zeros(Nl,1); pval_bis = zeros(Nl,1); pval_bsp = zeros(Nl,1);
    hval_bip = zeros(Nl,1); hval_bis = zeros(Nl,1); hval_bsp = zeros(Nl,1);

    med_i = zeros(Nl,NM);
    med_p = zeros(Nl,NM);
    med_s = zeros(Nl,NM);
    iqr_i = zeros(Nl,NM);
    iqr_p = zeros(Nl,NM);
    iqr_s = zeros(Nl,NM);

    best_i = zeros(Nl,J);
    best_s = zeros(Nl,J);
    best_p = zeros(Nl,J);

    N = 12e3;
    Signals = zeros(J,Nl,N);
    True = zeros(J,Nl,N);
    S_Spl = zeros(J,Nl,N);
    S_Pch = zeros(J,Nl,N);

    Err_Imp = zeros(Nl,J);
    Err_Spl = zeros(Nl,J);
    Err_Pch = zeros(Nl,J);
    for j=1:J
        name_r = files_r(j).name;
        load(fullfile(drt_r,name_r))

        err_imp(j,:,:) = St.Err_Imp(:,:,Ei)./max(St.True(1,:));
        err_spl(j,:,:) = St.Err_Spl(:,:,Ei)./max(St.True(1,:));
        err_pch(j,:,:) = St.Err_Pch(:,:,Ei)./max(St.True(1,:));

    end

    for i=1:Nl

        med_i(i,:) = median(squeeze(err_imp(:,i,:)),1,'omitnan');
        med_s(i,:) = median(squeeze(err_spl(:,i,:)),1,'omitnan');
        med_p(i,:) = median(squeeze(err_pch(:,i,:)),1,'omitnan');

        [~,best_i(i,:)] = min(squeeze(err_imp(:,i,:)),[],2);
        [~,best_s(i,:)] = min(squeeze(err_spl(:,i,:)),[],2);
        [~,best_p(i,:)] = min(squeeze(err_pch(:,i,:)),[],2);

        for j=1:J
            best_err_imp(j) = squeeze(err_imp(j,i,best_i(i,j)));
            best_err_spl(j) = squeeze(err_spl(j,i,best_i(i,j)));
            best_err_pch(j) = squeeze(err_pch(j,i,best_i(i,j)));
        end

        times_imp(i,:) = mean(St.Times_Imp,'omitnan');
        times_decomp(i,:) = mean(St.Times_Decomp,'omitnan');
        times_pch(i,:) = mean(St.Times_Pch,'omitnan');
        times_spl(i,:) = mean(St.Times_Spl,'omitnan');

        med_ib(i) = median(best_err_imp);
        med_sb(i) = median(best_err_spl);
        med_pb(i) = median(best_err_pch);

        [pval_bip(i),hval_bip(i)] = signrank(best_err_imp,best_err_pch,'alpha',alpha_b);
        [pval_bis(i),hval_bis(i)] = signrank(best_err_imp,best_err_spl,'alpha',alpha_b);
        [pval_bsp(i),hval_bsp(i)] = signrank(best_err_spl,best_err_pch,'alpha',alpha_b);

        times_imp(i,:) = mean(St.Times_Imp,'omitnan');
        times_decomp(i,:) = mean(St.Times_Decomp,'omitnan');
        times_pch(i,:) = mean(St.Times_Pch,'omitnan');
        times_spl(i,:) = mean(St.Times_Spl,'omitnan');

        Err_Imp(i,:) = best_err_imp;
        Err_Spl(i,:) = best_err_spl;
        Err_Pch(i,:) = best_err_pch;

    end
    PVAL_BIP(k,:) = pval_bip;
    PVAL_BIS(k,:) = pval_bis;
    PVAL_BSP(k,:) = pval_bsp;
    HVAL_BIP(k,:) = hval_bip;
    HVAL_BIS(k,:) = hval_bis;
    HVAL_BSP(k,:) = hval_bsp;
    MED_IB(k,:) = med_ib;
    MED_SB(k,:) = med_sb;
    MED_PB(k,:) = med_pb;


    h = {Err_Imp',Err_Spl',Err_Pch'};
    figure(1)
    set(gcf,'Position',[729 162 1192 800])
    subplot(2,3,k)
    set(gca,'Position',PosB(k,:))
    bx = boxplotGroup(h,'PrimaryLabels',{'BestImp','HaLI[BI](s)','HaLI[BI](p)'},'SecondaryLabels',{'5%','10%','15%','20%'});
    text(LabX(k),LabY(k),label,'FontSize',14,'interpreter','latex')
    ylabel('NMAE')
    ylim(YLim(k,:))
    %title(label)
    for aux=1:length(bx.axis2.XTickLabel)
        bx.axis2.XTickLabel{aux} = ['\newline\newline\newline' bx.axis2.XTickLabel{aux} ];
    end

    figure(2)
    set(gcf,'Position',[847 262 1074 700])
    subplot(2,3,k)
    catM = discretize(best_i,1:NM+1,'categorical',ImpNames);
    histogram(catM,'DisplayOrder','descend')
    text(LabX2(k),LabY2(k),label,'FontSize',14,'interpreter','latex')
    %histogram(catM)
    %title([label '. Best Imputation'])
end
