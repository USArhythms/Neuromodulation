% calculate subject averages

NE_order = order(NE_Idx);

tmp = struct;
tmp.Hb.LR_perf = f_regImages(Hb_model.Hb_LR_perf,refParcellation,settings.parcellation,0).*BM;
tmp.Hb.IRFx2_perf = f_regImages(Hb_model.Hb_IRFx2_perf,refParcellation,settings.parcellation,0).*BM;
tmp.HbO.LR_perf = f_regImages(Hb_model.HbO_LR_perf,refParcellation,settings.parcellation,0).*BM;
tmp.HbO.IRFx2_perf = f_regImages(Hb_model.HbO_IRFx2_perf,refParcellation,settings.parcellation,0).*BM;
tmp.Hb.LR_A = f_regImages(Hb_model.Hb_LR_A,refParcellation,settings.parcellation,0).*BM;
tmp.Hb.LR_B = f_regImages(Hb_model.Hb_LR_B,refParcellation,settings.parcellation,0).*BM;
tmp.Hb.IRFx2_A = f_regImages(Hb_model.Hb_IRFx2_A,refParcellation,settings.parcellation,0).*BM;
tmp.Hb.IRFx2_B = f_regImages(Hb_model.Hb_IRFx2_B,refParcellation,settings.parcellation,0).*BM;
tmp.HbO.LR_A = f_regImages(Hb_model.HbO_LR_A,refParcellation,settings.parcellation,0).*BM;
tmp.HbO.LR_B = f_regImages(Hb_model.HbO_LR_B,refParcellation,settings.parcellation,0).*BM;
tmp.HbO.IRFx2_A = f_regImages(Hb_model.HbO_IRFx2_A,refParcellation,settings.parcellation,0).*BM;
tmp.HbO.IRFx2_B = f_regImages(Hb_model.HbO_IRFx2_B,refParcellation,settings.parcellation,0).*BM;

subAvg.Hb.LR_perf = NaN(500,600,numel(NE_order));
subAvg.Hb.IRFx2_perf = NaN(500,600,numel(NE_order));
subAvg.Fig2.LR_perf = NaN(500,600,numel(NE_order));
subAvg.Fig2.IRFx2_perf = NaN(500,600,numel(NE_order));
subAvg.Hb.LR_A = NaN(500,600,numel(NE_order));
subAvg.Hb.LR_B = NaN(500,600,numel(NE_order));
subAvg.Hb.IRFx2_A = NaN(500,600,numel(NE_order));
subAvg.Hb.IRFx2_B = NaN(500,600,numel(NE_order));
subAvg.Hb.tA = NaN(numel(NE_order),1);
subAvg.Hb.tB = NaN(numel(NE_order),1);
subAvg.Hb.IRFx2_IRF = NaN(151,2,numel(NE_order));
subAvg.HbO.LR_A = NaN(500,600,numel(NE_order));
subAvg.HbO.LR_B = NaN(500,600,numel(NE_order));
subAvg.HbO.IRFx2_A = NaN(500,600,numel(NE_order));
subAvg.HbO.IRFx2_B = NaN(500,600,numel(NE_order));
subAvg.HbO.tA = NaN(numel(NE_order),1);
subAvg.HbO.tB = NaN(numel(NE_order),1);
subAvg.HbO.IRFx2_IRF = NaN(151,2,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.Hb.LR_perf(:,:,i) = mean(tmp.Hb.LR_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Hb.IRFx2_perf(:,:,i) = mean(tmp.Hb.IRFx2_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.HbO.LR_perf(:,:,i) = mean(tmp.HbO.LR_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.HbO.IRFx2_perf(:,:,i) = mean(tmp.HbO.IRFx2_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Hb.LR_A(:,:,i) = mean(tmp.Hb.LR_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Hb.LR_B(:,:,i) = mean(tmp.Hb.LR_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Hb.IRFx2_A(:,:,i) = mean(tmp.Hb.IRFx2_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Hb.IRFx2_B(:,:,i) = mean(tmp.Hb.IRFx2_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Hb.LR_tA(i) = mean([Hb_model.Hb_LR_tA{NE_order(i).Runs}]);
    subAvg.Hb.LR_tB(i) = mean([Hb_model.Hb_LR_tB{NE_order(i).Runs}]);
    subAvg.Hb.IRFx2_IRF(:,:,i) = mean(cat(3,Hb_model.Hb_IRFx2_IRF{NE_order(i).Runs}),3);
    subAvg.HbO.LR_A(:,:,i) = mean(tmp.HbO.LR_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.HbO.LR_B(:,:,i) = mean(tmp.HbO.LR_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.HbO.IRFx2_A(:,:,i) = mean(tmp.HbO.IRFx2_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.HbO.IRFx2_B(:,:,i) = mean(tmp.HbO.IRFx2_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.HbO.LR_tA(i) = mean([Hb_model.HbO_LR_tA{NE_order(i).Runs}]);
    subAvg.HbO.LR_tB(i) = mean([Hb_model.HbO_LR_tB{NE_order(i).Runs}]);
    subAvg.HbO.IRFx2_IRF(:,:,i) = mean(cat(3,Hb_model.HbO_IRFx2_IRF{NE_order(i).Runs}),3);
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

%% Fig HbO A

f = figure;
f_plotMap(mean(subAvg.HbO.LR_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR A',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

f = figure;
f_plotMap(mean(subAvg.HbO.LR_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR B',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig HbO B

barData = {};
barData{1} = subAvg.HbO.LR_tA;
barData{2} = subAvg.HbO.LR_tB;

f = figure;
f_plotBar(barData,colors=[c_Ca;c_GRAB],legend={'tA','tB'},ylabel='r',title='Timing Coefficients');

%% Fig HbO C

f = figure;
f_plotMap(mean(subAvg.HbO.LR_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='LR Performance',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig HbO D

f = figure;
f_plotMap((mean(subAvg.HbO.LR_perf,3,'omitnan')-mean(subAvg.Fig2.g_LR_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR vs. Global IRF',clabel='\Deltar');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig HbR E

f = figure;
f_plotMap(mean(subAvg.Hb.LR_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR A',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

f = figure;
f_plotMap(mean(subAvg.Hb.LR_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR B',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig HbR F

barData = {};
barData{1} = subAvg.Hb.LR_tA;
barData{2} = subAvg.Hb.LR_tB;

f = figure;
f_plotBar(barData,colors=[c_Ca;c_GRAB],legend={'tA','tB'},ylabel='r',title='Timing Coefficients');

%% Fig HbR G

f = figure;
f_plotMap(mean(subAvg.Hb.LR_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='LR Performance',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig HbR H

f = figure;
f_plotMap((mean(subAvg.Hb.LR_perf,3,'omitnan')-mean(subAvg.Fig2.g_LR_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR vs. Global IRF',clabel='\Deltar');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig HbO I

f = figure;
f_plotMap(mean(subAvg.HbO.IRFx2_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 A',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

f = figure;
f_plotMap(mean(subAvg.HbO.IRFx2_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 B',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig HbO J

f = figure;
meanSig = mean(subAvg.HbO.IRFx2_IRF(:,1,:),3);
error = std(subAvg.HbO.IRFx2_IRF(:,1,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig,error,color=c_Ca);
meanSig = mean(subAvg.HbO.IRFx2_IRF(:,2,:),3);
error = std(subAvg.HbO.IRFx2_IRF(:,2,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig,error,color=c_GRAB);
xlim([-2 7]);
xlabel('Time (s)');
ylabel('a.u.');
legend('','IRF(t_0^A,\tau_A)','','IRF(t_0^B,\tau_B)');
title('Double IRF');
set(gca,'FontSize',14);

%% Fig HbR K

f = figure;
f_plotMap(mean(subAvg.HbO.IRFx2_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='IRFx2 Performance',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig HbO L

f = figure;
f_plotMap((mean(subAvg.HbO.IRFx2_perf,3,'omitnan')-mean(subAvg.Fig2.g_IRFx2_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 vs. Global IRF',clabel='\Deltar');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig HbR M

f = figure;
f_plotMap(mean(subAvg.Hb.IRFx2_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 A',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

f = figure;
f_plotMap(mean(subAvg.Hb.IRFx2_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 B',clabel='');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig HbR N

f = figure;
meanSig = mean(subAvg.Hb.IRFx2_IRF(:,1,:),3);
error = std(subAvg.Hb.IRFx2_IRF(:,1,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig,error,color=c_Ca);
meanSig = mean(subAvg.Hb.IRFx2_IRF(:,2,:),3);
error = std(subAvg.Hb.IRFx2_IRF(:,2,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig,error,color=c_GRAB);
xlim([-2 7]);
xlabel('Time (s)');
ylabel('a.u.');
legend('','IRF(t_0^A,\tau_A)','','IRF(t_0^B,\tau_B)');
title('Double IRF');
set(gca,'FontSize',14);

%% Fig HbR O

f = figure;
f_plotMap(mean(subAvg.Hb.IRFx2_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='IRFx2 Performance',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig HbR P

f = figure;
f_plotMap((mean(subAvg.Hb.IRFx2_perf,3,'omitnan')-mean(subAvg.Fig2.g_IRFx2_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 vs. Global IRF',clabel='\Deltar');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%% Fig HbO Q

barData = {};
barData{1} = squeeze(mean(subAvg.Fig1.inv_perf.*plotBM,[1,2],'omitnan'));
barData{2} = squeeze(mean(subAvg.Fig1.SSp_perf.*plotBM,[1,2],'omitnan'));
barData{3} = squeeze(mean(subAvg.Fig1.var_perf.*plotBM,[1,2],'omitnan'));
barData{4} = squeeze(mean(subAvg.Fig2.g_LR_perf.*plotBM,[1,2],'omitnan'));
barData{5} = squeeze(mean(subAvg.Fig2.g_IRFx2_perf.*plotBM,[1,2],'omitnan'));
barData{6} = squeeze(mean(subAvg.HbO.LR_perf.*plotBM,[1,2],'omitnan'));
barData{7} = squeeze(mean(subAvg.HbO.IRFx2_perf.*plotBM,[1,2],'omitnan'));
barData{8} = squeeze(mean(subAvg.Hb.LR_perf.*plotBM,[1,2],'omitnan'));
barData{9} = squeeze(mean(subAvg.Hb.IRFx2_perf.*plotBM,[1,2],'omitnan'));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=[repmat(c_Yellow,3,1);repmat(c_darkCyan,2,1);repmat(c_Ca,2,1);repmat(c_pupil,2,1)],legend={'Invariant','SSp','Variant','LR','IRFx2','HbO LR','HbO IRFx2','HbR LR','HbR IRFx2'},ylabel='r',title='Model Performance Comparison');

[h,p] = f_kstest(barData,0.01);