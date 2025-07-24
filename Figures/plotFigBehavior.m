% calculate subject averages

NE_order = order(NE_Idx);

tmp.Beh.R_rfp_HD_HbT = f_regImages(Behavior.R.rfp_HD_low_HbT_low,refParcellation,settings.parcellation,0).*BM;
tmp.Beh.R_gfp_HD_HbT = f_regImages(Behavior.R.gfp_HD_low_HbT_low,refParcellation,settings.parcellation,0).*BM;
tmp.Beh.R_rfp_HD_gfp_HD = f_regImages(Behavior.R.rfp_HD_low_gfp_HD_low,refParcellation,settings.parcellation,0).*BM;
tmp.Beh.NE_IRF_perf = f_regImages(Behavior.NE_IRF_perf,refParcellation,settings.parcellation,0).*BM;
tmp.Beh.SPG_rfp_HD = f_adjust_SPG(Behavior.SPG.rfp_HD,4097,1);
tmp.Beh.SPG_gfp_HD = f_adjust_SPG(Behavior.SPG.gfp_HD,4097,1);
tmp.Beh.SPG_HbT = f_adjust_SPG(Behavior.SPG.HbT,4097,1);
tmp.Beh.COH_rfp_HD_gfp_HD = f_adjust_SPG(Behavior.COH.rfp_HD_gfp_HD,4097,0);
tmp.Beh.COH_rfp_HD_HbT = f_adjust_SPG(Behavior.COH.rfp_HD_HbT,4097,0);
tmp.Beh.COH_gfp_HD_HbT = f_adjust_SPG(Behavior.COH.gfp_HD_HbT,4097,0);
tmp.Beh.PHI_rfp_HD_gfp_HD = f_adjust_SPG(Behavior.PHI.rfp_HD_gfp_HD,4097,0);
tmp.Beh.PHI_rfp_HD_HbT = f_adjust_SPG(Behavior.PHI.rfp_HD_HbT,4097,0);
tmp.Beh.PHI_gfp_HD_HbT = f_adjust_SPG(Behavior.PHI.gfp_HD_HbT,4097,0);

subAvg.Beh.XC_gfp_HD_HbT = NaN(201,M);
subAvg.Beh.XC_rfp_HD_HbT = NaN(201,M);
subAvg.Beh.XC_rfp_HD_gfp_HD = NaN(201,M);
subAvg.Beh.R_rfp_HD_HbT = NaN(500,600,M);
subAvg.Beh.R_gfp_HD_HbT = NaN(500,600,M);
subAvg.Beh.R_rfp_HD_gfp_HD = NaN(500,600,M);
subAvg.Beh.R_behavior = NaN(6,6,M);
subAvg.Beh.SPG_rfp_HD = NaN(4097,M);
subAvg.Beh.SPG_gfp_HD = NaN(4097,M);
subAvg.Beh.SPG_HbT = NaN(4097,M);
subAvg.Beh.COH_rfp_HD_gfp_HD = NaN(4097,M);
subAvg.Beh.COH_rfp_HD_HbT = NaN(4097,M);
subAvg.Beh.COH_gfp_HD_HbT = NaN(4097,M);
subAvg.Beh.PHI_rfp_HD_gfp_HD = NaN(4097,M);
subAvg.Beh.PHI_rfp_HD_HbT = NaN(4097,M);
subAvg.Beh.PHI_gfp_HD_HbT = NaN(4097,M);
subAvg.Beh.NE_IRF_perf = NaN(500,600,M);
subAvg.Beh.NE_IRF = NaN(151,M);
subAvg.Beh.GRAB_conn = NaN(12,12,M);
subAvg.Beh.XC_gfp_HD_HbT_allen = NaN(201,12,M);

for i = 1:M
    subAvg.Beh.XC_gfp_HD_HbT(:,i) = mean(cat(2,Behavior.XC.gfp_HD_HbT{order(i).Runs}),2);
    subAvg.Beh.XC_rfp_HD_HbT(:,i) = mean(cat(2,Behavior.XC.rfp_HD_HbT{order(i).Runs}),2);
    subAvg.Beh.XC_rfp_HD_gfp_HD(:,i) = mean(cat(2,Behavior.XC.rfp_HD_gfp_HD{order(i).Runs}),2);
    subAvg.Beh.R_rfp_HD_HbT(:,:,i) = mean(tmp.Beh.R_rfp_HD_HbT(:,:,order(i).Runs),3,'omitnan');
    subAvg.Beh.R_gfp_HD_HbT(:,:,i) = mean(tmp.Beh.R_gfp_HD_HbT(:,:,order(i).Runs),3,'omitnan');
    subAvg.Beh.R_rfp_HD_gfp_HD(:,:,i) = mean(tmp.Beh.R_rfp_HD_gfp_HD(:,:,order(i).Runs),3,'omitnan');
    subAvg.Beh.R_behavior(:,:,i) = mean(cat(3,Behavior.R.signals{order(i).Runs}),3,'omitnan');
    subAvg.Beh.SPG_rfp_HD(:,i) = mean(cat(2,tmp.Beh.SPG_rfp_HD{order(i).Runs}),2);
    subAvg.Beh.SPG_gfp_HD(:,i) = mean(cat(2,tmp.Beh.SPG_gfp_HD{order(i).Runs}),2);
    subAvg.Beh.SPG_HbT(:,i) = mean(cat(2,tmp.Beh.SPG_HbT{order(i).Runs}),2);
    subAvg.Beh.COH_rfp_HD_gfp_HD(:,i) = mean(cat(2,tmp.Beh.COH_rfp_HD_gfp_HD{order(i).Runs}),2);
    subAvg.Beh.COH_rfp_HD_HbT(:,i) = mean(cat(2,tmp.Beh.COH_rfp_HD_HbT{order(i).Runs}),2);
    subAvg.Beh.COH_gfp_HD_HbT(:,i) = mean(cat(2,tmp.Beh.COH_gfp_HD_HbT{order(i).Runs}),2);
    subAvg.Beh.PHI_rfp_HD_gfp_HD(:,i) = mean(cat(2,tmp.Beh.PHI_rfp_HD_gfp_HD{order(i).Runs}),2);
    subAvg.Beh.PHI_rfp_HD_HbT(:,i) = mean(cat(2,tmp.Beh.PHI_rfp_HD_HbT{order(i).Runs}),2);
    subAvg.Beh.PHI_gfp_HD_HbT(:,i) = mean(cat(2,tmp.Beh.PHI_gfp_HD_HbT{order(i).Runs}),2);
    subAvg.Beh.NE_IRF_perf(:,:,i) = mean(tmp.Beh.NE_IRF_perf(:,:,order(i).Runs),3,'omitnan');
    try 
        subAvg.Beh.NE_IRF(:,i) = mean(cat(2,Behavior.NE_IRF_IRF{order(i).Runs}),2);
        subAvg.Beh.XC_gfp_HD_HbT_allen(:,:,i) = mean(cat(3,GRAB_FC.gfp_HD_vs_HbT_low{order(i).Runs}),3);
    end
    subAvg.Beh.GRAB_conn(:,:,i) = mean(cat(3,GRAB_FC.FC_detrend{order(i).Runs}),3);
end

fr = Behavior.SPG.fr;

plotBM = refBM;
plotBM(:,1:300) = NaN;

%% Fig Behavior A

f = figure;
meanSig = mean(subAvg.Beh.SPG_rfp_HD,2);
error = std(subAvg.Beh.SPG_rfp_HD,0,2)/sqrt(M);
f_plotLineError(fr,meanSig,error,color=c_Ca,log=1,lineWidth=3);
meanSig = mean(subAvg.Beh.SPG_gfp_HD(:,ACh_Idx),2);
error = std(subAvg.Beh.SPG_gfp_HD(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(fr,meanSig,error,color=c_Orange,log=1,lineWidth=3);
meanSig = mean(subAvg.Beh.SPG_gfp_HD(:,NE_Idx),2);
error = std(subAvg.Beh.SPG_gfp_HD(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(fr,meanSig,error,color=c_GRAB,log=1,lineWidth=3);
meanSig = mean(subAvg.Beh.SPG_HbT,2);
error = std(subAvg.Beh.SPG_HbT,0,2)/sqrt(M);
f_plotLineError(fr,meanSig,error,color=c_HbT,log=1,lineWidth=3);

xlim([0.05,5]);
xlabel('F (Hz)');
ylabel('Normalized Power');
legend('','Ca^2^+','','ACh','','NE','','HbT');
set(gca,'FontSize',14);

%% Fig Behavior B

t = 10:-0.1:-10;

f = figure;
meanSig = mean(subAvg.Beh.XC_rfp_HD_gfp_HD(:,NE_Idx),2);
error = std(subAvg.Beh.XC_rfp_HD_gfp_HD(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(t,meanSig,error,color=c_GRAB,lineWidth=3);
meanSig = mean(subAvg.Beh.XC_rfp_HD_gfp_HD(:,ACh_Idx),2);
error = std(subAvg.Beh.XC_rfp_HD_gfp_HD(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(t,meanSig,error,color=c_Orange,lineWidth=3);
xlim([-5 5]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);
title('x vs. Ca^2^+');
legend('','NE','','ACh');

%% Fig Behavior C

R_Beh = mean(subAvg.Beh.R_behavior,3);
gfp_R = mean(subAvg.Beh.R_behavior(:,:,ACh_Idx),3);
R_Beh(2,:) = gfp_R(2,:);
R_Beh(:,2) = gfp_R(2,:);
gfp_R = mean(subAvg.Beh.R_behavior(:,:,NE_Idx),3);
R_Beh(4:7,:) = R_Beh(3:6,:);
R_Beh(:,4:7) = R_Beh(:,3:6);
R_Beh(3,[1,3,4,5,6,7]) = gfp_R(2,:);
R_Beh([1,3,4,5,6,7],3) = gfp_R(2,:);
R_Beh(2,3) = 0;R_Beh(3,2) = 0;

barData = {};
barData{1} = squeeze(subAvg.Beh.R_behavior(2,4,ACh_Idx));
barData{2} = squeeze(subAvg.Beh.R_behavior(2,5,ACh_Idx));
barData{3} = squeeze(subAvg.Beh.R_behavior(2,6,ACh_Idx));
barData{4} = squeeze(subAvg.Beh.R_behavior(2,4,NE_Idx));
barData{5} = squeeze(subAvg.Beh.R_behavior(2,5,NE_Idx));
barData{6} = squeeze(subAvg.Beh.R_behavior(2,6,NE_Idx));
barData{7} = squeeze(subAvg.Beh.R_behavior(1,4,:));
barData{8} = squeeze(subAvg.Beh.R_behavior(1,5,:));
barData{9} = squeeze(subAvg.Beh.R_behavior(1,6,:));
barData{10} = squeeze(subAvg.Beh.R_behavior(3,4,:));
barData{11} = squeeze(subAvg.Beh.R_behavior(3,5,:));
barData{12} = squeeze(subAvg.Beh.R_behavior(3,6,:));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=repmat([c_pupil;0,0.7,0.7;0,0,0],4,1),legend={'Pupil diameter','Whisking','Movement'},ylabel='r',title='Model Performance Comparison',ylim=[0,0.8]);

f = figure;
imagesc(R_Beh);
colormap cmpbbr;
clim([-1 1]);
axis image off;
c = colorbar;
c.Label.String = 'r';
set(gca,'FontSize',14);

%% Fig Behavior D

f = figure;
f_plotMap(mean(subAvg.Beh.R_rfp_HD_gfp_HD(:,:,ACh_Idx),3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='ACh vs. Ca^2^+',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

barData = {};
barData{1} = squeeze(mean(subAvg.Beh.R_rfp_HD_gfp_HD(:,:,ACh_Idx).*plotBM,[1,2],'omitnan'));

f = figure(Position=[100,100,200,330]);
[dataMean, dataSEM] = f_plotBar(barData,colors=c_darkCyan,ylabel='r');
ylim([0 1]);

%% Fig Behavior E

f = figure(Position=[100,100,450,450]);
meanSig = mean(subAvg.Beh.COH_rfp_HD_gfp_HD(:,ACh_Idx),2);
error = std(subAvg.Beh.COH_rfp_HD_gfp_HD(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(fr,meanSig,error,color=c_Orange,log=1,lineWidth=3);
ylim([0,1]);
set(gca,'YScale','linear','FontSize',14);
ylabel('Coherence');

yyaxis right;
meanSig = mean(subAvg.Beh.PHI_rfp_HD_gfp_HD(:,ACh_Idx),2);
error = std(subAvg.Beh.PHI_rfp_HD_gfp_HD(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(fr,meanSig,error,color=[0,0,0],log=1,lineWidth=3);
ylim(pi*[-1 1]);
xlim([0.05,5]);
set(gca,'YScale','linear','FontSize',14);
xlabel('F (Hz)');
ylabel('Phi (rad)');

%% Fig Behavior G

mIdx = 3;
runIdx = order(mIdx).Runs;

cmp = cmpinf;
cmp = cmp(1:end-20,:);

X = Fig1.GRAB_global(runIdx);
Y = Fig1.SSp_perf_dt(runIdx);
f = figure;
f_multiScatter(X,Y,cmp=cmp,alpha=0.5,lineWidth=2,xlim=[-6 8],ylim=[-1 1],ylabel='r',xlabel='ACh');

barData = {};
barData{1} = squeeze(subAvg.Fig1.SSp_perf_vs_GRAB_global(2,ACh_Idx));

f = figure(Position=[100,100,200,330]);
[dataMean, dataSEM] = f_plotBar(barData,colors=c_darkCyan,ylabel='r');
ylim([-1 1]);

%% Fig Behavior H

f = figure;
f_plotMap(mean(subAvg.Beh.R_rfp_HD_gfp_HD(:,:,NE_Idx),3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='NE vs. Ca^2^+',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

barData = {};
barData{1} = squeeze(mean(subAvg.Beh.R_rfp_HD_gfp_HD(:,:,NE_Idx).*plotBM,[1,2],'omitnan'));

f = figure(Position=[100,100,200,330]);
[dataMean, dataSEM] = f_plotBar(barData,colors=c_darkCyan,ylabel='r');
ylim([0 1]);

%% Fig Behavior I

f = figure(Position=[100,100,450,450]);
meanSig = mean(subAvg.Beh.COH_rfp_HD_gfp_HD(:,NE_Idx),2);
error = std(subAvg.Beh.COH_rfp_HD_gfp_HD(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(fr,meanSig,error,color=c_GRAB,log=1,lineWidth=3);
ylim([0,1]);
set(gca,'YScale','linear','FontSize',14);
ylabel('Coherence');

yyaxis right;
meanSig = mean(subAvg.Beh.PHI_rfp_HD_gfp_HD(:,NE_Idx),2);
error = std(subAvg.Beh.PHI_rfp_HD_gfp_HD(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(fr,meanSig,error,color=[0,0,0],log=1,lineWidth=3);
ylim(pi*[-1 1]);
xlim([0.05,5]);
set(gca,'YScale','linear','FontSize',14);
xlabel('F (Hz)');
ylabel('Phi (rad)');

%% Fig Behavior J

meanSig = mean(subAvg.Beh.NE_IRF(:,NE_Idx),2);
error = std(subAvg.Beh.NE_IRF(:,NE_Idx),0,2)/sqrt(M_NE);

f = figure(Position=[100,100,300,200]);
f_plotLineError(-5:0.1:10,meanSig,error,color=c_darkCyan);
xlim([-3 7]);
xlabel('Time (s)');
ylabel('a.u.');
set(gca,'FontSize',14);
title('NE IRF');

f = figure;
f_plotMap(mean(subAvg.Beh.NE_IRF_perf(:,:,NE_Idx),3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='NE IRF vs. Ca^2^+',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

barData = {};
barData{1} = squeeze(mean(subAvg.Beh.NE_IRF_perf(:,:,NE_Idx).*plotBM,[1,2],'omitnan'));

f = figure(Position=[100,100,200,330]);
[dataMean, dataSEM] = f_plotBar(barData,colors=c_darkCyan,ylabel='r');
ylim([0 1]);

%% Fig Behavior K

f = figure;
f_plotFC(mean(subAvg.Beh.GRAB_conn(:,:,NE_Idx),3),1,cmp=cmpvir,bounds=[0 1],title='NE Connectivity',clabel='r');

%% Fig Behavior L

f = figure(Position=[100,100,450,450]);
meanSig = mean(subAvg.Beh.COH_rfp_HD_HbT,2);
error = std(subAvg.Beh.COH_rfp_HD_HbT,0,2)/sqrt(M);
f_plotLineError(fr,meanSig,error,color=c_Ca,log=1,lineWidth=3);
ylim([0,1]);
set(gca,'YScale','linear','FontSize',14);
ylabel('Coherence');

yyaxis right;
meanSig = mean(subAvg.Beh.PHI_rfp_HD_HbT,2);
error = std(subAvg.Beh.PHI_rfp_HD_HbT,0,2)/sqrt(M);
f_plotLineError(fr,meanSig,error,color=[0,0,0],log=1,lineWidth=3);
ylim(pi*[-1 1]);
xlim([0.05,5]);
set(gca,'YScale','linear','FontSize',14);
xlabel('F (Hz)');
ylabel('Phi (rad)');

f = figure(Position=[100,100,450,450]);
meanSig = mean(subAvg.Beh.COH_gfp_HD_HbT(:,NE_Idx),2);
error = std(subAvg.Beh.COH_gfp_HD_HbT(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(fr,meanSig,error,color=c_GRAB,log=1,lineWidth=3);
ylim([0,1]);
set(gca,'YScale','linear','FontSize',14);
ylabel('Coherence');

yyaxis right;
meanSig = mean(subAvg.Beh.PHI_gfp_HD_HbT(:,NE_Idx),2);
error = std(subAvg.Beh.PHI_gfp_HD_HbT(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(fr,meanSig,error,color=[0,0,0],log=1,lineWidth=3);
ylim(pi*[-1 1]);
xlim([0.05,5]);
set(gca,'YScale','linear','FontSize',14);
xlabel('F (Hz)');
ylabel('Phi (rad)');

f = figure(Position=[100,100,450,450]);
meanSig = mean(subAvg.Beh.COH_gfp_HD_HbT(:,ACh_Idx),2);
error = std(subAvg.Beh.COH_gfp_HD_HbT(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(fr,meanSig,error,color=c_Orange,log=1,lineWidth=3);
ylim([0,1]);
set(gca,'YScale','linear','FontSize',14);
ylabel('Coherence');

yyaxis right;
meanSig = mean(subAvg.Beh.PHI_gfp_HD_HbT(:,ACh_Idx),2);
error = std(subAvg.Beh.PHI_gfp_HD_HbT(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(fr,meanSig,error,color=[0,0,0],log=1,lineWidth=3);
ylim(pi*[-1 1]);
xlim([0.05,5]);
set(gca,'YScale','linear','FontSize',14);
xlabel('F (Hz)');
ylabel('Phi (rad)');

%% Fig Behavior M

t = 10:-0.1:-10;

f = figure;
meanSig = mean(subAvg.Beh.XC_gfp_HD_HbT(:,NE_Idx),2);
error = std(subAvg.Beh.XC_gfp_HD_HbT(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(t,meanSig,error,color=c_GRAB);
meanSig = mean(subAvg.Beh.XC_gfp_HD_HbT(:,ACh_Idx),2);
error = std(subAvg.Beh.XC_gfp_HD_HbT(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(t,meanSig,error,color=c_Orange);
meanSig = mean(subAvg.Beh.XC_rfp_HD_HbT,2);
error = std(subAvg.Beh.XC_rfp_HD_HbT,0,2)/sqrt(M);
f_plotLineError(t,meanSig,error,color=c_Ca);
xlim([-5 5]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);
title('x vs. HbT');
legend('','NE','','ACh','','Ca^2^+');

%% plot 2P xcorr

t = 10:-0.1:-10;

f = figure;
meanSig = mean(subAvg.Beh.XC_gfp_HD_HbT_allen(:,2,NE_Idx),3);
error = std(subAvg.Beh.XC_gfp_HD_HbT_allen(:,2,NE_Idx),0,3)/sqrt(M_NE);
f_plotLineError(t,meanSig,error,color=c_GRAB,lineWidth=3);
xlim([-5 5]);
ylim(0.7*[-1 1]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);
title('NE(MOs) vs. HbT');

f = figure;
meanSig = mean(subAvg.Beh.XC_gfp_HD_HbT_allen(:,5,NE_Idx),3);
error = std(subAvg.Beh.XC_gfp_HD_HbT_allen(:,5,NE_Idx),0,3)/sqrt(M_NE);
f_plotLineError(t,meanSig,error,color=c_GRAB,lineWidth=3);
xlim([-5 5]);
ylim(0.7*[-1 1]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);
title('NE(SSp-ll) vs. HbT');

f = figure;
meanSig = mean(subAvg.Beh.XC_gfp_HD_HbT_allen(:,12,NE_Idx),3);
error = std(subAvg.Beh.XC_gfp_HD_HbT_allen(:,12,NE_Idx),0,3)/sqrt(M_NE);
f_plotLineError(t,meanSig,error,color=c_GRAB,lineWidth=3);
xlim([-5 5]);
ylim(0.7*[-1 1]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);
title('NE(VISp) vs. HbT');

%%
function spg = f_adjust_SPG(data,length,norm)
    N = numel(data);
    spg = data;
    for i = 1:N
        if numel(spg{i}) ~= length
            spg{i} = movmean(spg{i},2);
            spg{i} = spg{i}(1:2:end);
        end
        if norm
            spg{i} = spg{i}/mean(spg{i})/5;
        end
    end
end