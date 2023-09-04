addpath(genpath('C:\Users\Intel\Dropbox\Github'))
drt = '/media/Datos/joaquinruiz/MissingDatSynth/MultiInt_Noise2';

N = 4000;
fs = 4000;
load(fullfile(drt,['Results_MissingDataMultiInt_200_10.mat']))
ratio = 0.05:0.05:0.2;
Nl = length(ratio);
vSNR = [10, 20, Inf];
ErrorCrit = upper(cellstr(S.Metrics));
Ei = 1;
ImpNames = cellstr(S.ImpNames);
NM = length(ImpNames);
J = 100;

%%
med_imp = zeros(Nl,NM);med_spl = med_imp;med_pch = med_imp;med_lin = med_imp;
best_imp = zeros(Nl,J);best_spl = best_imp;best_pch = best_imp; best_lin = best_imp;
pval_ip = zeros(Nl,NM);pval_is = pval_ip;pval_il = pval_ip;pval_sp = pval_ip;pval_sl = pval_ip;pval_pl = pval_ip;
hval_ip = zeros(Nl,NM);hval_is = hval_ip;hval_il = hval_ip;hval_sp = hval_ip;hval_sl = hval_ip;hval_pl = hval_ip;
times_imp = zeros(Nl,NM);times_decomp = times_imp; times_pch = times_imp;
times_spl = times_imp; times_lin = times_imp;
err_imp_best = zeros(J,1); err_spl_best = err_imp_best; err_pch_best = err_imp_best;

med_ib = zeros(Nl,1); med_sb = med_ib; med_pb = med_ib;
pval_bip = zeros(Nl,1); pval_bis = zeros(Nl,1); pval_bsp = zeros(Nl,1);
hval_bip = zeros(Nl,1); hval_bis = zeros(Nl,1); hval_bsp = zeros(Nl,1);
alpha_b = 0.05/3;

YLim(:,:,1) = [-0.013 0.875; -0.07 1.57; -0.07 1.38; -0.12 1.40];
YLim(:,:,2) = [-0.08 0.64; -0.08 1.08; -0.08 1.23; -0.08 1.4];
BP = [0.03,0.03,0,0; 0 0 0 0]';
UP = [0.8, 1.43 ,1.25, 1.27; 0.58 0.98 1.13 1.27]';
PosBox = [0.11 0.61 0.34 0.33; 0.55 0.61 0.34 0.33; 0.11 0.18 0.34 0.33; 0.55 0.18 0.34 0.33];
edgec = [0.4 0.4 0.4];
for k=1:length(vSNR)
    SNR = vSNR(k);
    for i=1:Nl
        L = round(ratio(i)*N);
        load(fullfile(drt,['Results_MissingDataMultiInt_' num2str(L) '_' num2str(SNR) '.mat']))

        err_imp = S.Err_Imp(:,:,Ei);
        err_spl = S.Err_Spl(:,:,Ei);
        err_pch = S.Err_Pch(:,:,Ei);
        err_lin = S.Err_Lin(:,:,Ei);

        err = {err_imp,err_spl,err_pch,err_lin};
        %figure(i)
        %boxplotGroup(err,'PrimaryLabels',{'Imp','Spl','Pch','Lin'},'SecondaryLabels',ImpNames)
        %title(ErroCrit{Ei})
        %hold off


        for j=1:NM
            [pval_ip(i,j),hval_ip(i,j)] = signrank(err_imp(:,j),err_pch(:,j),'alpha',alpha_b);
            [pval_is(i,j),hval_is(i,j)] = signrank(err_imp(:,j),err_spl(:,j),'alpha',alpha_b);
            %[pval_il(i,j),hval_il(i,j)] = signrank(err_imp(:,j),err_lin(:,j),'alpha',alpha_b);
            [pval_sp(i,j),hval_sp(i,j)] = signrank(err_spl(:,j),err_pch(:,j),'alpha',alpha_b);
            %[pval_sl(i,j),hval_sl(i,j)] = signrank(err_spl(:,j),err_lin(:,j),'alpha',alpha_b);
            %[pval_pl(i,j),hval_pl(i,j)] = signrank(err_pch(:,j),err_lin(:,j),'alpha',alpha_b);

        end
        med_imp(i,:) = median(err_imp,1,'omitnan');
        med_spl(i,:) = median(err_spl,1,'omitnan');
        med_pch(i,:) = median(err_pch,1,'omitnan');
        med_lin(i,:) = median(err_lin,1,'omitnan');

        times_imp(i,:) = mean(S.Times_Imp,'omitnan');
        times_decomp(i,:) = mean(S.Times_Decomp,'omitnan');
        times_pch(i,:) = mean(S.Times_Pch,'omitnan');
        times_spl(i,:) = mean(S.Times_Spl,'omitnan');
        times_lin(i,:) = mean(S.Times_Lin,'omitnan');

        [~,best_imp(i,:)] = min(err_imp,[],2);
        [~,best_spl(i,:)] = min(err_spl,[],2);
        [~,best_pch(i,:)] = min(err_pch,[],2);
        [~,best_lin(i,:)] = min(err_lin,[],2);

        figure(k)
        subplot(2,2,i)
        catM = discretize(best_imp(i,:),[1:NM+1],'categorical',ImpNames);
        histogram(catM,'DisplayOrder','descend')
        xlabel([num2str(ratio(i)*100) '% Missing Data'],'FontSize',12)
        %title('Best Initial Imputation')

        for j=1:J
            err_imp_best(j) = err_imp(j,best_imp(i,j));
            err_spl_best(j) = err_spl(j,best_imp(i,j));
            err_pch_best(j) = err_pch(j,best_imp(i,j));
        end

        med_ib(i) = median(err_imp_best,'omitnan');
        med_sb(i) = median(err_spl_best,'omitnan');
        med_pb(i) = median(err_pch_best,'omitnan');


        [pval_bip(i),hval_bip(i)] = signrank(err_imp_best,err_pch_best,'alpha',alpha_b);
        [pval_bis(i),hval_bis(i)] = signrank(err_imp_best,err_spl_best,'alpha',alpha_b);
        [pval_bsp(i),hval_bsp(i)] = signrank(err_spl_best,err_pch_best,'alpha',alpha_b);
        Err_Imp_Best(:,i) = err_imp_best;
        Err_Spl_Best(:,i) = err_spl_best;
        Err_Pch_Best(:,i) = err_pch_best;
        figure(k)
        set(gcf,'Position',[902 284 1019 662])
        hax(i) = subplot(2,2,i);
        boxplot([err_imp_best err_spl_best err_pch_best])
        xticklabels({'BestImp','HaLI[BI](s)','HaLI[BI](p)'})
        xlabel([num2str(ratio(i)*100) '% Missing Data'],'FontSize',12)
        hax(i).TickLabelInterpreter = 'tex';
        set(hax(i),'Position',PosBox(i,:))

        ylabel(ErrorCrit{Ei})
        %text(1,BP(i,k),'|','Color',edgec)
        %text(1.95,BP(i,k),'|','Color',edgec)
        %line([1.01,1.96],[BP(i,k),BP(i,k)],'Color',edgec)
        %text(1.45,BP(i,k)+0.02*(YLim(i,2,k)-YLim(i,1,k)),'*','FontSize',16)

        %text(2.01,BP(i,k),'|','Color',edgec)
        %text(2.96,BP(i,k),'|','Color',edgec)
        %line([2.03,2.97],[BP(i,k),BP(i,k)],'Color',edgec)
        %text(2.45,BP(i,k)+0.02*(YLim(i,2,k)-YLim(i,1,k)),'*','FontSize',16)

        %text(0.99,UP(i,k),'|','Color',edgec)
        %text(2.99,UP(i,k),'|','Color',edgec)
        %line([1.01,3],[UP(i,k),UP(i,k)],'Color',edgec)
        %text(1.96,UP(i,k)+0.02*(YLim(i,2,k)-YLim(i,1,k)),'*','FontSize',16)

        %ylim(YLim(i,:,k))
    end
end

%%
drt_nl  = '/media/Datos/joaquinruiz/MissingDatSynth/MultiInt_Simu18';
ratio = 0.2;
L = round(ratio*N);
load(fullfile(drt_nl,['Results_MissingDataMultiInt_' num2str(L) '_phival1.mat']))
i=1;
sc = S.S_MS(i,:);

s1_ms = sc;

fmax = 0.5;
sigma = 5e-5;
F1ms = STFT_Gauss(s1_ms,length(s1_ms),sigma,fmax);

t = 0:1/fs:(N-1)/fs;
f = 0:fs/N:fs*fmax-fs/N;

figure(3)
set(gcf,'Position',[603 101 1026 821])
subplot(6,1,1)
plot(t,s1_ms)
set(gca,'FontSize',12)
%title('Missing Rate 20%')
text(0,3.5,'Noiseless','FontSize',14,'Interpreter','latex')
set(gca,'Position',[0.13 0.84 0.775 0.0964])
subplot(6,1,2:3)
%spectro(F1ms,t,f)
imagesc(abs(F1ms)), axis xy, colormap(1-gray)
set(gca,'FontSize',12)
%xlabel('Time [s]')
ylabel('Frequency [Hz]')
ylim([0 400])
set(gca,'Position',[0.13 0.553 0.775 0.2451])
    
%%
ratio = 0.2;
SNR = 10;
L = round(ratio*N);
load(fullfile(drt,['Results_MissingDataMultiInt_' num2str(L) '_' num2str(SNR) '.mat']))
i=20;


sc = S.S_MS(i,:);
s2_ms = sc;

fmax = 0.5;
sigma = 5e-5;
F2ms = STFT_Gauss(s2_ms,length(s2_ms),sigma,fmax);

t = 0:1/fs:(N-1)/fs;
f = 0:fs/N:fs*fmax-fs/N;

figure(3)
subplot(6,1,4)
plot(t,s2_ms)
set(gca,'FontSize',12)
yticks([-4 0 4])
%title('Missing Rate 20%')
text(0,5.7,'SNR = 10 dB','FontSize',14,'Interpreter','latex')
set(gca,'Position',[0.13 0.375 0.775 0.0928])
subplot(6,1,5:6)
spectro(F2ms,t,f)
imagesc(abs(F2ms)), axis xy, colormap(1-gray)
set(gca,'FontSize',12)
xlabel('Time [s]')
ylabel('Frequency [Hz]')
ylim([0 400])
set(gca,'Position',[0.13 0.08 0.775 0.2541])
