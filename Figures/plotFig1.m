
% calculate subject averages

tmp = struct;
tmp.inv_perf = f_regImages(Fig1.inv_perf,refParcellation,settings.parcellation,0).*BM;
tmp.SSp_perf = f_regImages(Fig1.SSp_perf,refParcellation,settings.parcellation,0).*BM;
tmp.var_perf = f_regImages(Fig1.var_perf,refParcellation,settings.parcellation,0).*BM;
tmp.inv_IRF = cat(2,Fig1.inv_IRF{:});
tmp.SSp_IRF = cat(2,Fig1.SSp_IRF{:});
tmp.var_IRF = cat(3,Fig1.var_IRF{:});
tmp.SSp_perf_vs_GRAB = cat(2,Fig1.SSp_perf_vs_GRAB{:});

subAvg.Fig1.inv_perf = NaN(500,600,M);
subAvg.Fig1.SSp_perf = NaN(500,600,M);
subAvg.Fig1.var_perf = NaN(500,600,M);
subAvg.Fig1.inv_IRF = NaN(101,M);
subAvg.Fig1.SSp_IRF = NaN(101,M);
subAvg.Fig1.var_IRF = NaN(101,12,M);
subAvg.Fig1.SSp_perf_vs_GRAB = NaN(12,M);
subAvg.Fig1.SSp_perf_vs_GRAB_global = NaN(12,M);

for i = 1:M
    subAvg.Fig1.inv_perf(:,:,i) = mean(tmp.inv_perf(:,:,order(i).Runs),3,'omitnan');
    subAvg.Fig1.SSp_perf(:,:,i) = mean(tmp.SSp_perf(:,:,order(i).Runs),3,'omitnan');
    subAvg.Fig1.var_perf(:,:,i) = mean(tmp.var_perf(:,:,order(i).Runs),3,'omitnan');
    subAvg.Fig1.inv_IRF(:,i) = mean(tmp.inv_IRF(:,order(i).Runs),2,'omitnan');
    subAvg.Fig1.SSp_IRF(:,i) = mean(tmp.SSp_IRF(:,order(i).Runs),2,'omitnan');
    subAvg.Fig1.var_IRF(:,:,i) = mean(tmp.var_IRF(:,:,order(i).Runs),3,'omitnan');
    subAvg.Fig1.SSp_perf_vs_GRAB(:,i) = mean(tmp.SSp_perf_vs_GRAB(:,order(i).Runs),2,'omitnan');
    subAvg.Fig1.SSp_perf_vs_GRAB_global(:,i) = mean(cat(1,Fig1.SSp_perf_vs_GRAB_global{order(i).Runs}),1,'omitnan');
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

%% Fig 1C

f = figure;
f_plotMap(mean(subAvg.Fig1.inv_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='Global IRF',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[1 0 0]);
end

meanSig = mean(subAvg.Fig1.inv_IRF,2);
error = std(subAvg.Fig1.inv_IRF,0,2)/sqrt(M);

f = figure(Position=[100,100,300,200]);
f_plotLineError(0:0.1:10,meanSig,error,color=c_darkCyan);
xlim([0 7]);
xlabel('Time (s)');
ylabel('a.u.');
set(gca,'FontSize',14);
title('Global IRF');

%% Fig 1D

f = figure;
f_plotMap(mean(subAvg.Fig1.SSp_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='Global IRF',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end
f_plotAllenRegion(4,2,linewidth=3,color=[1 0 0]);
f_plotAllenRegion(5,2,linewidth=3,color=[1 0 0]);

meanSig = mean(subAvg.Fig1.SSp_IRF,2);
error = std(subAvg.Fig1.SSp_IRF,0,2)/sqrt(M);

f = figure(Position=[100,100,300,200]);
f_plotLineError(0:0.1:10,meanSig,error,color=c_darkCyan);
xlim([0 7]);
xlabel('Time (s)');
ylabel('a.u.');
set(gca,'FontSize',14);
title('SSp IRF');

%% Fig 1E

f = figure;
f_plotMap(mean(subAvg.Fig1.var_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='Global IRF',clabel='r');
for i = 1:12        
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end
f_plotAllenRegion(2,2,linewidth=3,color=c_darkCyan);
f_plotAllenRegion(4,2,linewidth=3,color=c_Orange);
f_plotAllenRegion(5,2,linewidth=3,color=c_Orange);

f = figure(Position=[100,100,300,200]);
meanSig = mean(subAvg.Fig1.var_IRF(:,2),3);
error = std(subAvg.Fig1.var_IRF(:,2,:),0,3)/sqrt(M);
f_plotLineError(0:0.1:10,meanSig,error,color=c_darkCyan);
meanSig = mean(subAvg.Fig1.var_IRF(:,[4,5],:),2);
error = std(meanSig,0,3)/sqrt(M);
meanSig = mean(meanSig,3);
f_plotLineError(0:0.1:10,meanSig,error,color=c_Orange);
xlim([0 7]);
xlabel('Time (s)');
ylabel('a.u.');
legend('','MOs','','SSp-tr/ll');
title('Variant IRF');
set(gca,'FontSize',14);

%% plot 1F

mIdx = 15;
runIdx = order(mIdx).Runs;

cmp = cmpinf;
cmp = cmp(1:end-20,:);

X = Fig1.GRAB_global(runIdx);
Y = Fig1.SSp_perf_dt(runIdx);
X(30) = [];
Y(30) = [];
f = figure;
f_multiScatter(X,Y,cmp=cmp,alpha=0.5,lineWidth=2,xlim=[-6 8],ylim=[-1 1],ylabel='r',xlabel='NE');

barData = {};
barData{1} = squeeze(subAvg.Fig1.SSp_perf_vs_GRAB_global(2,ACh_Idx));

%%
f = figure(Position=[100,100,200,330]);
[dataMean, dataSEM] = f_plotBar(barData,colors=c_darkCyan,ylabel='r');
ylim(0.5*[-1 1]);



%% plot 1G

f = figure;
f_plotAllenMap(mean(subAvg.Fig1.SSp_perf_vs_GRAB_global(:,NE_Idx),2),cmp=cmpbbr,cLabel='r',mask=plotBM,cRange=[-0.3, 0.3]);

%% plot 1H

f = figure;
f_plotAllenMap(std(subAvg.Fig1.SSp_perf_vs_GRAB_global(:,NE_Idx),0,2)/sqrt(sum(NE_Idx)),cmp=cmpinf,cLabel='SEM',mask=plotBM,cRange=[0, 0.1]);

%% plot 1I

barData = {};
barData{1} = squeeze(mean(subAvg.Fig1.inv_perf.*plotBM,[1,2],'omitnan'));
barData{2} = squeeze(mean(subAvg.Fig1.SSp_perf.*plotBM,[1,2],'omitnan'));
barData{3} = squeeze(mean(subAvg.Fig1.var_perf.*plotBM,[1,2],'omitnan'));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=c_darkCyan,legend={'Invariant','SSp','Variant'},ylabel='r',title='Model Performance Comparison');

p = f_SRtest(barData,0.05);

%%

cmp = cmpinf;
% cmp = cmp(1:end-20,:);
cmp = c_GRAB;
X = cell(1);
Y = cell(1);
X{1} = grab(1:6000,5);
Y{1} = Ca(1:6000,5);
X{2} = grab(1:6000,5);
Y{2} = Ca(1:6000,5);
f = figure(Position=[100 100 500 450]);
f_multiScatter(X,Y,cmp=cmp,alpha=0.1,lineWidth=2,xlim=[-6 8],ylim=[-10 15],ylabel='r',xlabel='NE');
axis off;

