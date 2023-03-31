addpath(genpath('/home/sentey/Dropbox/Github'))
addpath(genpath('C:\Users\Intel\Dropbox\Github'))
addpath(genpath('/home/jruiz'))
severity = 'Normal';

Signals = {'Plethysmogram2','Airflow2','NasalPressure2','Thorax'};
Labels = {'PPG','AF','NP','THO'};
K = length(Signals);

ratio = 0.05:0.05:0.2;
Nl = length(ratio);
ratio_s = {'5%','10%','15%','20%'};
alpha_b = 0.05/3;

pvalmse_ip = zeros(Nl,K);
pvalmse_is = zeros(Nl,K);
pvalmse_sp = zeros(Nl,K);
pvalmse_il = zeros(Nl,K);
pvalmse_sl = zeros(Nl,K);
pvalmse_pl = zeros(Nl,K);
hvalmse_ip = zeros(Nl,K);
hvalmse_is = zeros(Nl,K);
hvalmse_sp = zeros(Nl,K);
hvalmse_il = zeros(Nl,K);
hvalmse_sl = zeros(Nl,K);
hvalmse_pl = zeros(Nl,K);

pvalmae_ip = zeros(Nl,K);
pvalmae_is = zeros(Nl,K);
pvalmae_sp = zeros(Nl,K);
pvalmae_il = zeros(Nl,K);
pvalmae_sl = zeros(Nl,K);
pvalmae_pl = zeros(Nl,K);
hvalmae_ip = zeros(Nl,K);
hvalmae_is = zeros(Nl,K);
hvalmae_sp = zeros(Nl,K);
hvalmae_il = zeros(Nl,K);
hvalmae_sl = zeros(Nl,K);
hvalmae_pl = zeros(Nl,K);

medmse_i = zeros(Nl,K);
medmse_p = zeros(Nl,K);
medmse_s = zeros(Nl,K);
medmse_l = zeros(Nl,K);
iqrmse_i = zeros(Nl,K);
iqrmse_p = zeros(Nl,K);
iqrmse_s = zeros(Nl,K);
iqrmse_l = zeros(Nl,K);

medmae_i = zeros(Nl,K);
medmae_p = zeros(Nl,K);
medmae_s = zeros(Nl,K);
medmae_l = zeros(Nl,K);
iqrmae_i = zeros(Nl,K);
iqrmae_p = zeros(Nl,K);
iqrmae_s = zeros(Nl,K);
iqrmae_l = zeros(Nl,K);

NM = 6;
for k=1:K
    signal = Signals{k};
    label = Labels{k};
    
    %drt_r = ['E:\MissingDataReal\' signal];
    drt_r = ['/media/Datos/joaquinruiz/MissingDataReal/' signal];
    files_r = dir(fullfile(drt_r,severity));
    
    files_r = files_r(3:end);
    
    J = length(files_r);
    TMSE_Best = zeros(Nl,J);
    TMSE_Spl = zeros(Nl,J);
    TMSE_Pch = zeros(Nl,J);
    TMSE_Lin = zeros(Nl,J);
    TMAE_Best = zeros(Nl,J);
    TMAE_Spl = zeros(Nl,J);
    TMAE_Pch = zeros(Nl,J);
    TMAE_Lin = zeros(Nl,J);
    
    name_r = files_r(1).name;
    load(fullfile(drt_r,severity,name_r))
    
    N = size(St.True,2);
    S_IMP = zeros(J,Nl,N);
    True = zeros(J,Nl,N);
    S_Spl = zeros(J,Nl,N);
    S_Pch = zeros(J,Nl,N);
    BestM = zeros(Nl,J);

    for j=1:J
        name_r = files_r(j).name;
        load(fullfile(drt_r,severity,name_r))
        
        if size(St.True,1)==Nl
            True(j,:,:) = St.True;
            S_IMP(j,:,:) = St.S_Imp;
            S_Spl(j,:,:) = St.S_Spl;
            S_Pch(j,:,:) = St.S_Pch;
            TMSE_Best(:,j) = St.TMSE_Best/(max(St.True(:))-min(St.True(:)))^2;
            TMSE_Spl(:,j) = St.TMSE_Spl/(max(St.True(:))-min(St.True(:)))^2;
            TMSE_Pch(:,j) = St.TMSE_Pch/(max(St.True(:))-min(St.True(:)))^2;
            TMSE_Lin(:,j) = St.TMSE_Lin/(max(St.True(:))-min(St.True(:)));
            TMAE_Best(:,j) = St.TMAE_Best/(max(St.True(:))-min(St.True(:)));
            TMAE_Spl(:,j) = St.TMAE_Spl/(max(St.True(:))-min(St.True(:)));
            TMAE_Pch(:,j) = St.TMAE_Pch/(max(St.True(:))-min(St.True(:)));
            TMAE_Lin(:,j) = St.TMAE_Lin/(max(St.True(:))-min(St.True(:)));

            BestM(:,j) = St.BestM;
        else

            True(j,:,:) = St.True(1:2:7,:);
            S_IMP(j,:,:) = St.S_Imp(1:2:7,:);
            S_Spl(j,:,:) = St.S_Spl(1:2:7,:);
            S_Pch(j,:,:) = St.S_Pch(1:2:7,:);
            TMSE_Best(:,j) = St.TMSE_Best(1:2:7,:)/(max(St.True(:))-min(St.True(:)))^2;
            TMSE_Spl(:,j) =  St.TMSE_Spl(1:2:7,:)/(max(St.True(:))-min(St.True(:)))^2;
            TMSE_Pch(:,j) =  St.TMSE_Pch(1:2:7,:)/(max(St.True(:))-min(St.True(:)))^2;
            TMSE_Lin(:,j) =  St.TMSE_Lin(1:2:7,:)/(max(St.True(:))-min(St.True(:)));
            TMAE_Best(:,j) = St.TMAE_Best(1:2:7,:)/(max(St.True(:))-min(St.True(:)));
            TMAE_Spl(:,j) =  St.TMAE_Spl(1:2:7,:)/(max(St.True(:))-min(St.True(:)));
            TMAE_Pch(:,j) =  St.TMAE_Pch(1:2:7,:)/(max(St.True(:))-min(St.True(:)));
            TMAE_Lin(:,j) =  St.TMAE_Lin(1:2:7,:)/(max(St.True(:))-min(St.True(:)));
            BestM(:,j) = St.BestM(1:2:7);
        end
    end
    
    FrecM = zeros(6,Nl);
    for i=1:6
        FrecM(i,:) = sum(BestM==i,2);
    end

    for i=1:Nl
        [pvalmse_ip(i,k),hvalmse_ip(i,k)] = signrank(TMSE_Best(i,:),TMSE_Pch(i,:),'alpha',alpha_b);
        [pvalmse_is(i,k),hvalmse_is(i,k)] = signrank(TMSE_Best(i,:),TMSE_Spl(i,:),'alpha',alpha_b);
        [pvalmse_sp(i,k),hvalmse_sp(i,k)] = signrank(TMSE_Spl(i,:),TMSE_Pch(i,:),'alpha',alpha_b);
        [pvalmse_il(i,k),hvalmse_il(i,k)] = signrank(TMSE_Best(i,:),TMSE_Lin(i,:),'alpha',alpha_b);
        [pvalmse_sl(i,k),hvalmse_sl(i,k)] = signrank(TMSE_Spl(i,:),TMSE_Lin(i,:),'alpha',alpha_b);
        [pvalmse_pl(i,k),hvalmse_pl(i,k)] = signrank(TMSE_Pch(i,:),TMSE_Lin(i,:),'alpha',alpha_b);
        
        [pvalmae_ip(i,k),hvalmae_ip(i,k)] = signrank(TMAE_Best(i,:),TMAE_Pch(i,:),'alpha',alpha_b);
        [pvalmae_is(i,k),hvalmae_is(i,k)] = signrank(TMAE_Best(i,:),TMAE_Spl(i,:),'alpha',alpha_b);
        [pvalmae_sp(i,k),hvalmae_sp(i,k)] = signrank(TMAE_Spl(i,:),TMAE_Pch(i,:),'alpha',alpha_b);
        [pvalmae_il(i,k),hvalmae_il(i,k)] = signrank(TMAE_Best(i,:),TMAE_Lin(i,:),'alpha',alpha_b);
        [pvalmae_sl(i,k),hvalmae_sl(i,k)] = signrank(TMAE_Spl(i,:),TMAE_Lin(i,:),'alpha',alpha_b);
        [pvalmae_pl(i,k),hvalmae_pl(i,k)] = signrank(TMAE_Pch(i,:),TMAE_Lin(i,:),'alpha',alpha_b);
        
        medmse_i(i,k) = median(TMSE_Best(i,:));
        medmse_p(i,k) = median(TMSE_Pch(i,:));
        medmse_s(i,k) = median(TMSE_Spl(i,:));
        medmse_l(i,k) = median(TMSE_Lin(i,:));
        
        iqrmse_i(i,k) = iqr(TMSE_Best(i,:));
        iqrmse_p(i,k) = iqr(TMSE_Pch(i,:));
        iqrmse_s(i,k) = iqr(TMSE_Spl(i,:));
        iqrmse_l(i,k) = iqr(TMSE_Lin(i,:));
        
        medmae_i(i,k) = median(TMAE_Best(i,:));
        medmae_p(i,k) = median(TMAE_Pch(i,:));
        medmae_s(i,k) = median(TMAE_Spl(i,:));

        iqrmae_i(i,k) = iqr(TMAE_Best(i,:));
        iqrmae_p(i,k) = iqr(TMAE_Pch(i,:));
        iqrmae_s(i,k) = iqr(TMAE_Spl(i,:));

    end
    h = {TMSE_Best',TMSE_Spl',TMSE_Pch'};
    figure(1)
    subplot(2,2,k)
    boxplotGroup(h,'PrimaryLabel',{'Best I','S','P'},'SecondaryLabel',ratio_s)
    ylabel('NMSE')
    %xlabel('Missing rate')
    title(label)
    hold off

    %figure(1+k)
    subplot(2,2,k)
    bar(FrecM)
    legend('5%','10%','15%','20%')
    % legend('TLM','LSE','DMD','GPR','ARF','ARB','Orientation','vertical','Location','northeast')
    xlim([0.5067    7])
    %xticklabels({'5%','10%','15%','20%'})
    xticklabels({'TLM','LSE','DMD','GPR','ARF','ARB'})
    title(label)

    catM = discretize(BestM(:),NM,'categorical',{'TLM','LSE','DMD','GPR','ARF','ARB'});
    figure(6)
    subplot(2,2,k)
    histogram(catM,'DisplayOrder','descend')
    title(label,'FontSize',12)

    
end