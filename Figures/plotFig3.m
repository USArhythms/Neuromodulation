% calculate subject averages

NE_order = order(NE_Idx);

subAvg.Fig3.lowNE_FC_Ca = NaN(12,12,numel(NE_order));
subAvg.Fig3.lowNE_FC_HbT = NaN(12,12,numel(NE_order));
subAvg.Fig3.highNE_FC_Ca = NaN(12,12,numel(NE_order));
subAvg.Fig3.highNE_FC_HbT = NaN(12,12,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.Fig3.lowNE_FC_Ca(:,:,i) = mean(cat(3,Fig3.lowNE_Ca{NE_order(i).Runs}),3,'omitnan');
    subAvg.Fig3.lowNE_FC_HbT(:,:,i) = mean(cat(3,Fig3.lowNE_HbT{NE_order(i).Runs}),3,'omitnan');
    subAvg.Fig3.highNE_FC_Ca(:,:,i) = mean(cat(3,Fig3.highNE_Ca{NE_order(i).Runs}),3,'omitnan');
    subAvg.Fig3.highNE_FC_HbT(:,:,i) = mean(cat(3,Fig3.highNE_HbT{NE_order(i).Runs}),3,'omitnan');
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

%% Fig 3 G

f = figure;
f_plotFC(mean(subAvg.Fig3.lowNE_FC_Ca,3),1,cmp=cmpvir,bounds=[0 1],title='Low NE Ca++ Connectivity',clabel='r');

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_Ca,3),1,cmp=cmpvir,bounds=[0 1],title='High NE Ca++ Connectivity',clabel='r');

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_Ca,3)-mean(subAvg.Fig3.lowNE_FC_Ca,3),0,cmp=cmpbbr,bounds=0.2*[-1 1],title='High-Low NE',clabel='\Deltar');

%% Fig 3 H

f = figure;
f_plotFC(mean(subAvg.Fig3.lowNE_FC_HbT,3),1,cmp=cmpvir,bounds=[0 1],title='Low NE HbT Connectivity',clabel='r');

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_HbT,3),1,cmp=cmpvir,bounds=[0 1],title='High NE HbT Connectivity',clabel='r');

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_HbT,3)-mean(subAvg.Fig3.lowNE_FC_Ca,3),0,cmp=cmpbbr,bounds=0.2*[-1 1],title='High-Low NE',clabel='\Deltar');

%% Fig 3 I

allenCaLow = mean(subAvg.Fig3.lowNE_FC_Ca,3);
allenCaLow = allenCaLow([2,5,12],:);
allenCaHigh = mean(subAvg.Fig3.highNE_FC_Ca,3);
allenCaHigh = allenCaHigh([2,5,12],:);

for i = 1:3
    f = figure;f_plotAllenMap(allenCaLow(1,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);

end
