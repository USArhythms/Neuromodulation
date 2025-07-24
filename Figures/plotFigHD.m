

load('/projectnb/devorlab/bcraus/AnalysisCode/Neuromodulation/Figures/results/HD/HD_metadata.mat','metadata');

mice = fieldnames(metadata);

subAvg.HD.XC_gfp_HD_HbT_low = NaN(201,numel(mice));
subAvg.HD.R_gfp_HD_HbT_low = cell(3,1);

for i = 1:3
    subAvg.HD.XC_gfp_HD_HbT_low(:,i) = mean([metadata.(mice{i}).XC_gfp_HD_HbT_low],2);
    subAvg.HD.R_gfp_HD_HbT_low{i} = mean(cat(3,metadata.(mice{i}).R_gfp_HD_HbT_low),3);
end

%% plot xcorr

t = 10:-0.1:-10;

f = figure;
meanSig = mean(subAvg.HD.XC_gfp_HD_HbT_low,2);
error = std(subAvg.HD.XC_gfp_HD_HbT_low,0,2)/sqrt(numel(mice));
f_plotLineError(t,meanSig,error,color=0.5*[1,1,1],lineWidth=3);

%%
t = 10:-0.1:-10;

f = figure;
meanSig = mean(subAvg.Beh.XC_gfp_HD_HbT(:,NE_Idx),2);
error = std(subAvg.Beh.XC_gfp_HD_HbT(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(t,meanSig,error,color=c_GRAB);
meanSig = mean(subAvg.Beh.XC_gfp_HD_HbT(:,ACh_Idx),2);
error = std(subAvg.Beh.XC_gfp_HD_HbT(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(t,meanSig,error,color=c_Orange);
meanSig = mean(subAvg.HD.XC_gfp_HD_HbT_low,2);
error = std(subAvg.HD.XC_gfp_HD_HbT_low,0,2)/sqrt(3);
f_plotLineError(t,meanSig,error,color=0.5*[1,1,1],lineWidth=3);
xlim([-5 5]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);
title('x vs. HbT');
legend('','NE','','ACh','','Ca^2^+');

%%

f = figure;
f_plotMap(mean(subAvg.Beh.R_gfp_HD_HbT(:,:,NE_Idx),3,'omitnan').*plotBM,cmp=cmpbbr,bounds=1*[-1 1],title='NE vs. Ca^2^+',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%%