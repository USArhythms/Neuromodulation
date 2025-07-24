% calculate subject averages

NE_order = order(NE_Idx);

tmp = struct;
tmp.Unf.LR_perf = f_regImages(unfiltered.LR_perf,refParcellation,settings.parcellation,0).*BM;
tmp.Unf.IRFx2_perf = f_regImages(unfiltered.IRFx2_perf,refParcellation,settings.parcellation,0).*BM;
tmp.Shu.LR_perf = f_regImages(shuffled.LR_perf,refParcellation,settings.parcellation,0).*BM;
tmp.Shu.IRFx2_perf = f_regImages(shuffled.IRFx2_perf,refParcellation,settings.parcellation,0).*BM;
tmp.Unf.LR_A = f_regImages(unfiltered.LR_A,refParcellation,settings.parcellation,0).*BM;
tmp.Unf.LR_B = f_regImages(unfiltered.LR_B,refParcellation,settings.parcellation,0).*BM;
tmp.Unf.IRFx2_A = f_regImages(unfiltered.IRFx2_A,refParcellation,settings.parcellation,0).*BM;
tmp.Unf.IRFx2_B = f_regImages(unfiltered.IRFx2_B,refParcellation,settings.parcellation,0).*BM;
tmp.Shu.LR_A = f_regImages(shuffled.LR_A,refParcellation,settings.parcellation,0).*BM;
tmp.Shu.LR_B = f_regImages(shuffled.LR_B,refParcellation,settings.parcellation,0).*BM;
tmp.Shu.IRFx2_A = f_regImages(shuffled.IRFx2_A,refParcellation,settings.parcellation,0).*BM;
tmp.Shu.IRFx2_B = f_regImages(shuffled.IRFx2_B,refParcellation,settings.parcellation,0).*BM;

subAvg.Unf.LR_perf = NaN(500,600,numel(NE_order));
subAvg.Unf.IRFx2_perf = NaN(500,600,numel(NE_order));
subAvg.Fig2.LR_perf = NaN(500,600,numel(NE_order));
subAvg.Fig2.IRFx2_perf = NaN(500,600,numel(NE_order));
subAvg.Unf.LR_A = NaN(500,600,numel(NE_order));
subAvg.Unf.LR_B = NaN(500,600,numel(NE_order));
subAvg.Unf.IRFx2_A = NaN(500,600,numel(NE_order));
subAvg.Unf.IRFx2_B = NaN(500,600,numel(NE_order));
subAvg.Unf.tA = NaN(numel(NE_order),1);
subAvg.Unf.tB = NaN(numel(NE_order),1);
subAvg.Unf.IRFx2_IRF = NaN(151,2,numel(NE_order));
subAvg.Shu.LR_A = NaN(500,600,numel(NE_order));
subAvg.Shu.LR_B = NaN(500,600,numel(NE_order));
subAvg.Shu.IRFx2_A = NaN(500,600,numel(NE_order));
subAvg.Shu.IRFx2_B = NaN(500,600,numel(NE_order));
subAvg.Shu.tA = NaN(numel(NE_order),1);
subAvg.Shu.tB = NaN(numel(NE_order),1);
subAvg.Shu.IRFx2_IRF = NaN(151,2,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.Unf.LR_perf(:,:,i) = mean(tmp.Unf.LR_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Unf.IRFx2_perf(:,:,i) = mean(tmp.Unf.IRFx2_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Shu.LR_perf(:,:,i) = mean(tmp.Shu.LR_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Shu.IRFx2_perf(:,:,i) = mean(tmp.Shu.IRFx2_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Unf.LR_A(:,:,i) = mean(tmp.Unf.LR_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Unf.LR_B(:,:,i) = mean(tmp.Unf.LR_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Unf.IRFx2_A(:,:,i) = mean(tmp.Unf.IRFx2_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Unf.IRFx2_B(:,:,i) = mean(tmp.Unf.IRFx2_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Unf.LR_tA(i) = mean([unfiltered.LR_tA{NE_order(i).Runs}]);
    subAvg.Unf.LR_tB(i) = mean([unfiltered.LR_tB{NE_order(i).Runs}]);
    subAvg.Unf.IRFx2_IRF(:,:,i) = mean(cat(3,unfiltered.IRFx2_IRF{NE_order(i).Runs}),3);
    subAvg.Shu.LR_A(:,:,i) = mean(tmp.Shu.LR_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Shu.LR_B(:,:,i) = mean(tmp.Shu.LR_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Shu.IRFx2_A(:,:,i) = mean(tmp.Shu.IRFx2_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Shu.IRFx2_B(:,:,i) = mean(tmp.Shu.IRFx2_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Shu.LR_tA(i) = mean([shuffled.LR_tA{NE_order(i).Runs}]);
    subAvg.Shu.LR_tB(i) = mean([shuffled.LR_tB{NE_order(i).Runs}]);
    subAvg.Shu.IRFx2_IRF(:,:,i) = mean(cat(3,shuffled.IRFx2_IRF{NE_order(i).Runs}),3);
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

%% Fig unfiltered A

f = figure;
f_plotMap(mean(subAvg.Unf.LR_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR A',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

f = figure;
f_plotMap(mean(subAvg.Unf.LR_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR B',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig unfiltered B

barData = {};
barData{1} = subAvg.Unf.LR_tA;
barData{2} = subAvg.Unf.LR_tB;

f = figure;
f_plotBar(barData,colors=[c_Ca;c_GRAB],legend={'tA','tB'},ylabel='r',title='Timing Coefficients');

%% Fig unfiltered C

f = figure;
f_plotMap(mean(subAvg.Unf.LR_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='LR Performance',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig unfiltered D

f = figure;
f_plotMap((mean(subAvg.Unf.LR_perf,3,'omitnan')-mean(subAvg.Fig2.g_LR_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR vs. Global IRF',clabel='\Deltar');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig unfiltered E

f = figure;
f_plotMap(mean(subAvg.Unf.IRFx2_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 A',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

f = figure;
f_plotMap(mean(subAvg.Unf.IRFx2_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 B',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig unfiltered F

f = figure;
meanSig = mean(subAvg.Unf.IRFx2_IRF(:,1,:),3);
error = std(subAvg.Unf.IRFx2_IRF(:,1,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig,error,color=c_Ca);
meanSig = mean(subAvg.Unf.IRFx2_IRF(:,2,:),3);
error = std(subAvg.Unf.IRFx2_IRF(:,2,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig,error,color=c_GRAB);
xlim([-2 7]);
xlabel('Time (s)');
ylabel('a.u.');
legend('','IRF(t_0^A,\tau_A)','','IRF(t_0^B,\tau_B)');
title('Double IRF');
set(gca,'FontSize',14);

%% Fig unfiltered G

f = figure;
f_plotMap(mean(subAvg.Unf.IRFx2_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='IRFx2 Performance',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig unfiltered H

f = figure;
f_plotMap((mean(subAvg.Unf.IRFx2_perf,3,'omitnan')-mean(subAvg.Fig2.g_IRFx2_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 vs. Global IRF',clabel='\Deltar');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig shuffled I

f = figure;
f_plotMap(mean(subAvg.Shu.LR_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR A',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

f = figure;
f_plotMap(mean(subAvg.Shu.LR_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR B',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig shuffled J

barData = {};
barData{1} = subAvg.Shu.LR_tA;
barData{2} = subAvg.Shu.LR_tB;

f = figure;
f_plotBar(barData,colors=[c_Ca;c_GRAB],legend={'tA','tB'},ylabel='r',title='Timing Coefficients');

%% Fig shuffled K

f = figure;
f_plotMap((mean(subAvg.Shu.LR_perf,3,'omitnan')-mean(subAvg.Fig1.inv_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR vs. Global IRF',clabel='\Deltar');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig shuffled L

f = figure;
f_plotMap(mean(subAvg.Shu.IRFx2_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 A',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

f = figure;
f_plotMap(mean(subAvg.Shu.IRFx2_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 B',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig shuffled M

f = figure;
meanSig = mean(subAvg.Shu.IRFx2_IRF(:,1,:),3);
error = std(subAvg.Shu.IRFx2_IRF(:,1,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig,error,color=c_Ca);
meanSig = mean(subAvg.Shu.IRFx2_IRF(:,2,:),3);
error = std(subAvg.Shu.IRFx2_IRF(:,2,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig,error,color=c_GRAB);
xlim([-2 7]);
xlabel('Time (s)');
ylabel('a.u.');
legend('','IRF(t_0^A,\tau_A)','','IRF(t_0^B,\tau_B)');
title('Double IRF');
set(gca,'FontSize',14);

%% Fig shuffled N

f = figure;
f_plotMap((mean(subAvg.Shu.IRFx2_perf,3,'omitnan')-mean(subAvg.Fig1.inv_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 vs. Global IRF',clabel='\Deltar');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig shuffled O

barData = {};
barData{1} = squeeze(mean(subAvg.Fig1.inv_perf.*plotBM,[1,2],'omitnan'));
barData{2} = squeeze(mean(subAvg.Fig1.SSp_perf.*plotBM,[1,2],'omitnan'));
barData{3} = squeeze(mean(subAvg.Fig1.var_perf.*plotBM,[1,2],'omitnan'));
barData{4} = squeeze(mean(subAvg.Fig2.g_LR_perf.*plotBM,[1,2],'omitnan'));
barData{5} = squeeze(mean(subAvg.Fig2.g_IRFx2_perf.*plotBM,[1,2],'omitnan'));
barData{6} = squeeze(mean(subAvg.Unf.LR_perf.*plotBM,[1,2],'omitnan'));
barData{7} = squeeze(mean(subAvg.Unf.IRFx2_perf.*plotBM,[1,2],'omitnan'));
barData{8} = squeeze(mean(subAvg.Shu.LR_perf.*plotBM,[1,2],'omitnan'));
barData{9} = squeeze(mean(subAvg.Shu.IRFx2_perf.*plotBM,[1,2],'omitnan'));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=[repmat(c_Yellow,3,1);repmat(c_darkCyan,2,1);repmat(c_Ca,2,1);repmat(c_pupil,2,1)],legend={'Invariant','SSp','Variant','LR','IRFx2','unfiltered LR','unfiltered IRFx2','shuffled LR','shuffled IRFx2'},ylabel='r',title='Model Performance Comparison');

[h,p] = f_kstest(barData,0.01);