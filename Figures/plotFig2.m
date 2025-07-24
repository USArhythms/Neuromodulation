% calculate subject averages

NE_order = order(NE_Idx);

tmp = struct;
tmp.LR_perf = f_regImages(Fig2.LR_perf,refParcellation,settings.parcellation,0).*BM;
tmp.IRFx2_perf = f_regImages(Fig2.IRFx2_perf,refParcellation,settings.parcellation,0).*BM;
tmp.s_LR_perf = f_regImages(shuffled.LR_perf,refParcellation,settings.parcellation,0).*BM;
tmp.s_IRFx2_perf = f_regImages(shuffled.IRFx2_perf,refParcellation,settings.parcellation,0).*BM;
tmp.g_LR_perf = f_regImages(Fig2.global_LR_perf,refParcellation,settings.parcellation,0).*BM;
tmp.g_IRFx2_perf = f_regImages(Fig2.global_IRFx2_perf,refParcellation,settings.parcellation,0).*BM;
tmp.LR_A = f_regImages(Fig2.global_LR_A,refParcellation,settings.parcellation,0).*BM;
tmp.LR_B = f_regImages(Fig2.global_LR_B,refParcellation,settings.parcellation,0).*BM;
tmp.IRFx2_A = f_regImages(Fig2.global_IRFx2_A,refParcellation,settings.parcellation,0).*BM;
tmp.IRFx2_B = f_regImages(Fig2.global_IRFx2_B,refParcellation,settings.parcellation,0).*BM;

subAvg.Fig2.LR_perf = NaN(500,600,numel(NE_order));
subAvg.Fig2.IRFx2_perf = NaN(500,600,numel(NE_order));
subAvg.Fig2.s_LR_perf = NaN(500,600,numel(NE_order));
subAvg.Fig2.s_IRFx2_perf = NaN(500,600,numel(NE_order));
subAvg.Fig2.g_LR_perf = NaN(500,600,numel(NE_order));
subAvg.Fig2.g_IRFx2_perf = NaN(500,600,numel(NE_order));
subAvg.Fig2.LR_A = NaN(500,600,numel(NE_order));
subAvg.Fig2.LR_B = NaN(500,600,numel(NE_order));
subAvg.Fig2.IRFx2_A = NaN(500,600,numel(NE_order));
subAvg.Fig2.IRFx2_B = NaN(500,600,numel(NE_order));
subAvg.Fig2.tA = NaN(numel(NE_order),1);
subAvg.Fig2.tB = NaN(numel(NE_order),1);
subAvg.Fig2.IRFx2_IRF = NaN(151,2,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.Fig2.LR_perf(:,:,i) = mean(tmp.LR_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Fig2.IRFx2_perf(:,:,i) = mean(tmp.IRFx2_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Fig2.s_LR_perf(:,:,i) = mean(tmp.s_LR_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Fig2.s_IRFx2_perf(:,:,i) = mean(tmp.s_IRFx2_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Fig2.g_LR_perf(:,:,i) = mean(tmp.g_LR_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Fig2.g_IRFx2_perf(:,:,i) = mean(tmp.g_IRFx2_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Fig2.LR_A(:,:,i) = mean(tmp.LR_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Fig2.LR_B(:,:,i) = mean(tmp.LR_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Fig2.IRFx2_A(:,:,i) = mean(tmp.IRFx2_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Fig2.IRFx2_B(:,:,i) = mean(tmp.IRFx2_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Fig2.LR_tA(i) = mean([Fig2.global_LR_tA{NE_order(i).Runs}]);
    subAvg.Fig2.LR_tB(i) = mean([Fig2.global_LR_tB{NE_order(i).Runs}]);
    subAvg.Fig2.IRFx2_IRF(:,:,i) = mean(cat(3,Fig2.global_IRFx2_IRF{NE_order(i).Runs}),3);
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

%% Fig 2 B

f = figure;
f_plotMap(mean(subAvg.Fig2.LR_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR A',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

f = figure;
f_plotMap(mean(subAvg.Fig2.LR_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR B',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig 2 C

barData = {};
barData{1} = subAvg.Fig2.LR_tA;
barData{2} = subAvg.Fig2.LR_tB;

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=[c_Ca;c_GRAB],legend={'tA','tB'},ylabel='r',title='Timing Coefficients');

%% Fig 2 D

f = figure;
f_plotMap(mean(subAvg.Fig2.g_LR_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='LR Performance',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig 2 E

f = figure;
f_plotMap((mean(subAvg.Fig2.g_LR_perf,3,'omitnan')-mean(subAvg.Fig1.inv_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR vs. Global IRF',clabel='\Deltar');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig 2 F

f = figure;
f_plotMap(mean(subAvg.Fig2.IRFx2_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 A',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

f = figure;
f_plotMap(mean(subAvg.Fig2.IRFx2_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 B',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig 2 G

f = figure;
meanSig = mean(subAvg.Fig2.IRFx2_IRF(:,1,:),3);
error = std(subAvg.Fig2.IRFx2_IRF(:,1,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig,error,color=c_Ca);
meanSig = mean(subAvg.Fig2.IRFx2_IRF(:,2,:),3);
error = std(subAvg.Fig2.IRFx2_IRF(:,2,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig,error,color=c_GRAB);
xlim([-2 7]);
xlabel('Time (s)');
ylabel('a.u.');
legend('','IRF(t_0^A,\tau_A)','','IRF(t_0^B,\tau_B)');
title('Double IRF');
set(gca,'FontSize',14);

%% Fig 2 H

f = figure;
f_plotMap(mean(subAvg.Fig2.g_IRFx2_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='IRFx2 Performance',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig 2 I

f = figure;
f_plotMap((mean(subAvg.Fig2.g_IRFx2_perf,3,'omitnan')-mean(subAvg.Fig1.inv_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 vs. Global IRF',clabel='\Deltar');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig 2 J

f = figure;
f_plotMap(mean(subAvg.Fig2.s_LR_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='Shuffled LR Performance',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

f = figure;
f_plotMap(mean(subAvg.Fig2.s_IRFx2_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='Shuffled IRFx2 Performance',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig 2 K

barData = {};
barData{1} = squeeze(mean(subAvg.Fig1.inv_perf.*plotBM,[1,2],'omitnan'));
barData{2} = squeeze(mean(subAvg.Fig1.SSp_perf.*plotBM,[1,2],'omitnan'));
barData{3} = squeeze(mean(subAvg.Fig1.var_perf.*plotBM,[1,2],'omitnan'));
barData{4} = squeeze(mean(subAvg.Fig2.g_LR_perf.*plotBM,[1,2],'omitnan'));
barData{5} = squeeze(mean(subAvg.Fig2.g_IRFx2_perf.*plotBM,[1,2],'omitnan'));
barData{6} = squeeze(mean(subAvg.Fig2.s_LR_perf.*plotBM,[1,2],'omitnan'));
barData{7} = squeeze(mean(subAvg.Fig2.s_IRFx2_perf.*plotBM,[1,2],'omitnan'));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=[repmat(c_Yellow,3,1);repmat(c_darkCyan,2,1);repmat(c_pupil,2,1)],legend={'Invariant','SSp','Variant','LR','IRFx2','shuffled LR','shuffled IRFx2'},ylabel='r',title='Model Performance Comparison');

[h,p] = f_kstest(barData,0.01);