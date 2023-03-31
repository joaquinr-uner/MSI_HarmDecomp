drt = '/media/Datos/joaquinruiz/MissingDatSynth/MultiInt_Simu10';
%drt = '';
N = 8000;
ratio = 0.05:0.05:0.2;
Nl = length(ratio);

pval_ip = zeros(1,Nl);
pval_is = pval_ip;
pval_sp = pval_ip;
pval_il = pval_ip;
pval_sl = pval_ip;
pval_pl = pval_ip;
hval_ip = pval_ip;
hval_is = pval_ip;
hval_sp = pval_ip;
hval_il = pval_ip;
hval_sl = pval_ip;
hval_pl = pval_ip;
med_i = pval_ip;
med_p = pval_ip;
med_s = pval_ip;
med_l = pval_ip;
iqr_i = pval_ip;
iqr_p = pval_ip;
iqr_s = pval_ip;
iqr_l = pval_ip;

Matxy = [0.07,0.75,0.25,0.22;
         0.35,0.75,0.25,0.22;
         0.63,0.75,0.25,0.22;
         0.07,0.45,0.25,0.22;
         0.35,0.45,0.25,0.22;
         0.63,0.45,0.25,0.22;
         0.07,0.15,0.25,0.22];

Yastk = [-0.1,-0.3;
         -0.1,-0.3;
         -0.1,-0.3;
         -0.1,-0.3;
         -0.1,-0.3;
         -0.1,-0.3;
         -0.1,-0.3];

alpha_b = 0.05/3;
Nsbp = ceil(sqrt(Nl));

Nip = zeros(1,Nl);
Nis = zeros(1,Nl);
Nil = zeros(1,Nl);

for i=1:Nl
    L = round(ratio(i)*N);
    load(fullfile(drt,['Results_MissingDataMultiInt_' num2str(L) '_phival1.mat']))

    [pval_ip(i),hval_ip(i)] = signrank(S.TMSE_Best,S.TMSE_Pch,'alpha',alpha_b);
    [pval_is(i),hval_is(i)] = signrank(S.TMSE_Best,S.TMSE_Spl,'alpha',alpha_b); 
    [pval_sp(i),hval_sp(i)] = signrank(S.TMSE_Spl,S.TMSE_Pch,'alpha',alpha_b);
    [pval_il(i),hval_il(i)] = signrank(S.TMSE_Best,S.TMSE_Lin,'alpha',alpha_b);
    [pval_sl(i),hval_sl(i)] = signrank(S.TMSE_Spl,S.TMSE_Lin,'alpha',alpha_b);
    [pval_pl(i),hval_pl(i)] = signrank(S.TMSE_Pch,S.TMSE_Lin,'alpha',alpha_b);

    Nip(i) = sum(S.TMSE_Best>S.TMSE_Pch);
    Nis(i) = sum(S.TMSE_Best>S.TMSE_Spl);
    Nil(i) = sum(S.TMSE_Best>S.TMSE_Lin);
    med_i(i) = median(S.TMSE_Best);
    med_p(i) = median(S.TMSE_Pch);
    med_s(i) = median(S.TMSE_Spl);
    med_l(i) = median(S.TMSE_Lin);
    iqr_i(i) = iqr(S.TMSE_Best);
    iqr_p(i) = iqr(S.TMSE_Pch);
    iqr_s(i) = iqr(S.TMSE_Spl);
    iqr_l(i) = iqr(S.TMSE_Lin);
    NM = size(S.ImpNames,1);

    figure(1)
    subplot(Nsbp,Nsbp,i)
    boxplot([S.TMSE_Best,S.TMSE_Spl,S.TMSE_Pch])
    %title([Type ' Wave-Shape. L = ' num2str(L)], ['H_{ip} = ' num2str(round(hval_ip(i)))])
    xticklabels({'Best Imp','Spl','Pch'})
    xlabel([num2str(ratio(i)*100) '% Missing Data'],'FontSize',12)
    
    catM = discretize(S.BestM,NM,'categorical',{'TLM','LSE','DMD','GPR','ARF','ARB'});
    figure(2)
    subplot(Nsbp,Nsbp,i)
    histogram(catM,'DisplayOrder','descend')
    xlabel([num2str(ratio(i)*100) '% Missing Data'],'FontSize',12)
    %hold on
    %if pval_is<alpha_b
    %    plot([1,2],[-0.11,-0.11],'k')
    %    text(0.97,-0.1,'|','FontSize',12)
    %    text(1.45,-0.05,'*','FontSize',14)
    %    text(1.97,-0.1,'|','FontSize',12)
    %    set(gca,'Position',Matxy(i,:));
    %end
    % 
    %if pval_ip<alpha_b
    %    plot([1,3],[-0.35,-0.35],'k')
    %    text(0.97,-0.35,'|','FontSize',12)
    %    text(1.95,-0.3,'*','FontSize',14)
    %    text(2.97,-0.35,'|','FontSize',12)
    %    set(gca,'Position',Matxy(i,:));
    %end
    %yl = ylim;
    %ylim([-0.5,0.9355])
    %hold off
    %ylabel('MSE')

    %h = {S.MSE_Best,S.MSE_Spl,S.MSE_Pch};
    %figure(2)
    %subplot(3,3,i)
    %boxplotGroup(h,'PrimaryLabels',{'I','S','P'},'SecondaryLabels',{'I1','I2','I3'})
    %title([Type ' Wave-Shape. L = ' num2str(L)], ['H_{ip} = ' num2str(round(hval_ip(i)))])
    
    %figure(3)
    %subplot(3,3,i)
    %histogram(S.BestM)
    %xlabel(['L = ' num2str(L)])

    %figure(4)
    %subplot(3,3,i)
    %plot(S.True(1,:))
    %hold on
    %plot(S.S_Imp(1,:))
    %plot(S.S_Pch(1,:))

end

%% 0.05
ratio = 0.05;
L = round(ratio*N);
load(fullfile(drt,['Results_MissingDataMultiInt_' num2str(L) '_phival1.mat']))
[~,i] = max(S.TMSE_Best-S.TMSE_Pch);
sc = S.True(i,:);
s_imp = S.S_Imp(i,:);
s_pch = S.S_Pch(i,:);

r = (sc-s_imp~=0);

sti = find(diff(r)==1)+1;
edi = find(diff(r)==-1);

s_ms = sc;
s_ms(r==1) = 0;

fmax = 0.25;
sigma = 1e-4;
Fms = STFT_Gauss(s_ms,length(s_ms),sigma,fmax);

fs = 4000;
N = 8000;
t = 0:1/fs:(N-1)/fs;
f = 0:fs/N:fs*fmax-fs/N;

figure(5)
subplot(3,3,1)
plot(t,s_ms)
set(gca,'FontSize',12)
ylabel('Amplitude [A.U.]')
title('Missing Rate 5%')
subplot(3,3,[4,7])
spectro(Fms,t,f)
set(gca,'FontSize',12)
xlabel('Time [s]')
ylabel('Frequency [Hz]')

figure(6)
subplot(3,3,1)
plot(t(sti(1):edi(1)),sc(sti(1):edi(1)))
hold on
plot(t(sti(1):edi(1)),s_imp(sti(1):edi(1)),'*r')
plot(t(sti(1):edi(1)),s_pch(sti(1):edi(1)),'+k')
xlim([t(sti(1)) t(edi(1))])
hold off
subplot(3,3,2)
plot(t(sti(2):edi(2)),sc(sti(2):edi(2)))
hold on
plot(t(sti(2):edi(2)),s_imp(sti(2):edi(2)),'*r')
plot(t(sti(2):edi(2)),s_pch(sti(2):edi(2)),'+k')
xlim([t(sti(2)) t(edi(2))])
hold off
subplot(3,3,3)
plot(t(sti(3):edi(3)),sc(sti(3):edi(3)))
hold on
plot(t(sti(3):edi(3)),s_imp(sti(3):edi(3)),'*r')
plot(t(sti(3):edi(3)),s_pch(sti(3):edi(3)),'+k')
xlim([t(sti(3)) t(edi(3))])
hold off

%% 0.1
ratio = 0.1;
L = round(ratio*N);
load(fullfile(drt,['Results_MissingDataMultiInt_' num2str(L) '_phival1.mat']))
[~,i] = max(S.TMSE_Best-S.TMSE_Pch);
sc = S.True(i,:);
s_imp = S.S_Imp(i,:);
s_pch = S.S_Pch(i,:);

r = (sc-s_imp~=0);

sti = find(diff(r)==1)+1;
edi = find(diff(r)==-1);

s_ms = sc;
s_ms(r==1) = 0;

fmax = 0.25;
sigma = 1e-4;
Fms = STFT_Gauss(s_ms,length(s_ms),sigma,fmax);

fs = 4000;
N = 8000;
t = 0:1/fs:(N-1)/fs;
f = 0:fs/N:fs*fmax-fs/N;

figure(5)
subplot(3,3,2)
plot(t,s_ms)
set(gca,'FontSize',12)
title('Missing Rate 10%')
subplot(3,3,[5,8])
spectro(Fms,t,f)
set(gca,'FontSize',12)
xlabel('Time [s]')

figure(6)
subplot(3,3,4)
plot(t(sti(1):edi(1)),sc(sti(1):edi(1)))
hold on
plot(t(sti(1):edi(1)),s_imp(sti(1):edi(1)),'*r')
plot(t(sti(1):edi(1)),s_pch(sti(1):edi(1)),'+k')
xlim([t(sti(1)) t(edi(1))])
hold off
subplot(3,3,5)
plot(t(sti(2):edi(2)),sc(sti(2):edi(2)))
hold on
plot(t(sti(2):edi(2)),s_imp(sti(2):edi(2)),'*r')
plot(t(sti(2):edi(2)),s_pch(sti(2):edi(2)),'+k')
xlim([t(sti(2)) t(edi(2))])
hold off
subplot(3,3,6)
plot(t(sti(3):edi(3)),sc(sti(3):edi(3)))
hold on
plot(t(sti(3):edi(3)),s_imp(sti(3):edi(3)),'*r')
plot(t(sti(3):edi(3)),s_pch(sti(3):edi(3)),'+k')
xlim([t(sti(3)) t(edi(3))])
hold off

%% 0.2

ratio = 0.2;
L = round(ratio*N);
load(fullfile(drt,['Results_MissingDataMultiInt_' num2str(L) '_phival1.mat']))
[~,i] = max(S.TMSE_Best-S.TMSE_Pch);
sc = S.True(i,:);
s_imp = S.S_Imp(i,:);
s_pch = S.S_Pch(i,:);

r = sc-s_imp~=0;

sti = find(diff(r)==1)+1;
edi = find(diff(r)==-1);

s_ms = sc;
s_ms(r==1) = 0;

fmax = 0.25;
sigma = 1e-4;
Fms = STFT_Gauss(s_ms,length(s_ms),sigma,fmax);

fs = 4000;
N = 8000;
t = 0:1/fs:(N-1)/fs;
f = 0:fs/N:fs*fmax-fs/N;

figure(5)
subplot(3,3,3)
plot(t,s_ms)
set(gca,'FontSize',12)
title('Missing Rate 20%')
subplot(3,3,[6,9])
spectro(Fms,t,f)
set(gca,'FontSize',12)
xlabel('Time [s]')

figure(6)
subplot(3,3,7)
plot(t(sti(1):edi(1)),sc(sti(1):edi(1)))
hold on
plot(t(sti(1):edi(1)),s_imp(sti(1):edi(1)),'*r')
plot(t(sti(1):edi(1)),s_pch(sti(1):edi(1)),'+k')
xlim([t(sti(1)) t(edi(1))])
hold off
subplot(3,3,8)
plot(t(sti(2):edi(2)),sc(sti(2):edi(2)))
hold on
plot(t(sti(2):edi(2)),s_imp(sti(2):edi(2)),'*r')
plot(t(sti(2):edi(2)),s_pch(sti(2):edi(2)),'+k')
xlim([t(sti(2)) t(edi(2))])
hold off
subplot(3,3,9)
plot(t(sti(3):edi(3)),sc(sti(3):edi(3)))
hold on
plot(t(sti(3):edi(3)),s_imp(sti(3):edi(3)),'*r')
plot(t(sti(3):edi(3)),s_pch(sti(3):edi(3)),'+k')
xlim([t(sti(3)) t(edi(3))])
hold off

figure(5)
set(gcf,'Position',[4 45 1032 917])

