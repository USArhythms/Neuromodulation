% calculate subject averages

NE_order = order(NE_Idx);

subAvg.NE_reg.Ca_vs_HbT = NaN(12,12,numel(NE_order));
subAvg.NE_reg.Ca_vs_HbT_reg = NaN(12,12,numel(NE_order));
subAvg.NE_reg.global_Ca_vs_HbT_reg = NaN(12,12,numel(NE_order));
subAvg.NE_reg.lowNE_HbT_reg = NaN(12,12,numel(NE_order));
subAvg.NE_reg.highNE_HbT_reg = NaN(12,12,numel(NE_order));
subAvg.NE_reg.global_lowNE_HbT_reg = NaN(12,12,numel(NE_order));
subAvg.NE_reg.global_highNE_HbT_reg = NaN(12,12,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.NE_reg.Ca_vs_HbT(:,:,i) = mean(cat(3,NE_reg.FC_R_HbT{NE_order(i).Runs}),3);
    subAvg.NE_reg.Ca_vs_HbT_reg(:,:,i) = mean(cat(3,NE_reg.FC_R_HbT_reg{NE_order(i).Runs}),3);
    subAvg.NE_reg.global_Ca_vs_HbT_reg(:,:,i) = mean(cat(3,NE_reg.global_FC_R_HbT_reg{NE_order(i).Runs}),3);
    subAvg.NE_reg.lowNE_HbT_reg(:,:,i) = mean(cat(3,NE_reg.lowNE_FC{NE_order(i).Runs}),3);
    subAvg.NE_reg.highNE_HbT_reg(:,:,i) = mean(cat(3,NE_reg.highNE_FC{NE_order(i).Runs}),3);
    subAvg.NE_reg.global_lowNE_HbT_reg(:,:,i) = mean(cat(3,NE_reg.global_lowNE_FC{NE_order(i).Runs}),3);
    subAvg.NE_reg.global_highNE_HbT_reg(:,:,i) = mean(cat(3,NE_reg.global_highNE_FC{NE_order(i).Runs}),3);
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

%% Fig NE_reg B

f = figure;
f_plotFC(mean(subAvg.NE_reg.Ca_vs_HbT,3),0,cmp=cmpbbr,bounds=[-1 1],title='HbT',clabel='r');

f = figure;
f_plotFC(mean(subAvg.NE_reg.Ca_vs_HbT_reg,3),0,cmp=cmpbbr,bounds=[-1 1],title='NE \\ HbT',clabel='r');

f = figure;
f_plotFC(mean(subAvg.NE_reg.global_Ca_vs_HbT_reg,3),0,cmp=cmpbbr,bounds=[-1 1],title='Global NE \\ HbT',clabel='r');

%% Fig NE_reg C

f = figure;
f_plotFC(mean(subAvg.NE_reg.lowNE_HbT_reg,3),1,cmp=cmpvir,bounds=[0 1],title='Low NE Ca++ Connectivity',clabel='r');

f = figure;
f_plotFC(mean(subAvg.NE_reg.highNE_HbT_reg,3),1,cmp=cmpvir,bounds=[0 1],title='High NE Ca++ Connectivity',clabel='r');

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,NE_reg.lowNE_FC{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,NE_reg.highNE_FC{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

f = figure;
f_plotFC(mean(subAvg.NE_reg.highNE_HbT_reg,3)-mean(subAvg.NE_reg.lowNE_HbT_reg,3),0,cmp=cmpbbr,bounds=0.1*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));

%% Fig NE_reg D

f = figure;
f_plotFC(mean(subAvg.NE_reg.global_lowNE_HbT_reg,3),1,cmp=cmpvir,bounds=[0 1],title='Low NE Ca++ Connectivity',clabel='r');

f = figure;
f_plotFC(mean(subAvg.NE_reg.global_highNE_HbT_reg,3),1,cmp=cmpvir,bounds=[0 1],title='High NE Ca++ Connectivity',clabel='r');

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,NE_reg.global_lowNE_FC{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,NE_reg.global_highNE_FC{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

f = figure;
f_plotFC(mean(subAvg.NE_reg.global_highNE_HbT_reg,3)-mean(subAvg.NE_reg.global_lowNE_HbT_reg,3),0,cmp=cmpbbr,bounds=0.1*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));
