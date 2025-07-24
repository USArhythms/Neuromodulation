% calculate subject averages

NE_order = order(NE_Idx);

subAvg.Fig3.lowNE_FC_Ca = NaN(12,12,numel(NE_order));
subAvg.Fig3.lowNE_FC_HbT = NaN(12,12,numel(NE_order));
subAvg.Fig3.highNE_FC_Ca = NaN(12,12,numel(NE_order));
subAvg.Fig3.highNE_FC_HbT = NaN(12,12,numel(NE_order));
subAvg.Fig3.FC_Ca_HbT_vs_GRAB = NaN(numel(NE_order),1);
subAvg.Fig3.FC_Ca_vs_GRAB = NaN(12,12,numel(NE_order));
subAvg.Fig3.FC_HbT_vs_GRAB = NaN(12,12,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.Fig3.lowNE_FC_Ca(:,:,i) = mean(cat(3,Fig3.lowNE_Ca{NE_order(i).Runs}),3,'omitnan');
    subAvg.Fig3.lowNE_FC_HbT(:,:,i) = mean(cat(3,Fig3.lowNE_HbT{NE_order(i).Runs}),3,'omitnan');
    subAvg.Fig3.highNE_FC_Ca(:,:,i) = mean(cat(3,Fig3.highNE_Ca{NE_order(i).Runs}),3,'omitnan');
    subAvg.Fig3.highNE_FC_HbT(:,:,i) = mean(cat(3,Fig3.highNE_HbT{NE_order(i).Runs}),3,'omitnan');
    subAvg.Fig3.FC_Ca_HbT_vs_GRAB(i) = mean(cat(1,Fig3.FC_Ca_HbT_vs_GRAB{NE_order(i).Runs}));
    subAvg.Fig3.FC_Ca_vs_GRAB(:,:,i) = mean(cat(3,Fig3.FC_Ca_vs_GRAB{NE_order(i).Runs}),3);
    subAvg.Fig3.FC_HbT_vs_GRAB(:,:,i) = mean(cat(3,Fig3.FC_HbT_vs_GRAB{NE_order(i).Runs}),3);
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

%% Fig 3 D

mIdx = 8;
runIdx = order(mIdx).Runs;

cmp = cmpinf;
cmp = cmp(1:220,:);

X = Fig3.NE(runIdx);
Y = Fig3.FC_Ca_MOs_SSpll(runIdx);
f = figure;
f_multiScatter(X,Y,cmp=cmp,alpha=0.5,lineWidth=2,xlim=[-4 5],ylim=[0 1],ylabel='r',xlabel='NE');

Y = Fig3.FC_HbT_MOs_SSpll(runIdx);
f = figure;
f_multiScatter(X,Y,cmp=cmp,alpha=0.5,lineWidth=2,xlim=[-4 5],ylim=[0 1],ylabel='r',xlabel='NE');

X = cat(1,X{:});
f = figure;
h = histfit(X,30,'kernel');
h(1).FaceAlpha = 0.2;
h(1).FaceColor = [0 0 0];
h(1).EdgeColor = 'none';
h(2).Color = [0 0 0];
h(2).LineWidth = 1;
xlim([-4 5]);
axis off;

%% Fig 3 E
X = Fig3.NE(runIdx);
Y = Fig3.FC_Ca_vs_HbT(runIdx);
f = figure;
f_multiScatter(X,Y,cmp=cmp,alpha=0.5,lineWidth=2,xlim=[-4 5],ylim=[0 1],ylabel='r',xlabel='NE');

%% Fig 3 F

barData = {};
barData{1} = squeeze(subAvg.Fig3.FC_Ca_vs_GRAB(2,5,:));
barData{2} = squeeze(subAvg.Fig3.FC_HbT_vs_GRAB(2,5,:));
barData{3} = subAvg.Fig3.FC_Ca_HbT_vs_GRAB;

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=[c_Ca;c_HbT;c_darkCyan],legend={'Ca^2^+','HbT','Ca^2^+ vs. HbT'},ylabel='r');

%% Fig 3 G

f = figure;
f_plotFC(mean(subAvg.Fig3.lowNE_FC_Ca,3),1,cmp=cmpvir,bounds=[0 1],title='Low NE Ca++ Connectivity',clabel='r');

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_Ca,3),1,cmp=cmpvir,bounds=[0 1],title='High NE Ca++ Connectivity',clabel='r');

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,Fig3.lowNE_Ca{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,Fig3.highNE_Ca{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_Ca,3)-mean(subAvg.Fig3.lowNE_FC_Ca,3),0,cmp=cmpbbr,bounds=0.25*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));

%% Fig 3 H

f = figure;
f_plotFC(mean(subAvg.Fig3.lowNE_FC_HbT,3),1,cmp=cmpvir,bounds=[0 1],title='Low NE HbT Connectivity',clabel='r');

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_HbT,3),1,cmp=cmpvir,bounds=[0 1],title='High NE HbT Connectivity',clabel='r');

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,Fig3.lowNE_HbT{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,Fig3.highNE_HbT{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_HbT,3)-mean(subAvg.Fig3.lowNE_FC_HbT,3),0,cmp=cmpbbr,bounds=0.1*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));

%% Fig 3 I

regions = [2,5,12];

allenCaLow = mean(subAvg.Fig3.lowNE_FC_Ca,3);
allenCaLow = allenCaLow(regions,:);
allenCaHigh = mean(subAvg.Fig3.highNE_FC_Ca,3);
allenCaHigh = allenCaHigh(regions,:);

for i = 1:3
    f = figure;f_plotAllenMap(allenCaLow(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,sprintf('Fig3I_%01i.png',(i-1)*3+1),'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenCaHigh(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,sprintf('Fig3I_%01i.png',(i-1)*3+2),'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenCaHigh(i,:)-allenCaLow(i,:),cmp=cmpbbr,mask=plotBM,cRange=0.25*[-1,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=0.7*[1 1 1]);
    colorbar off;exportgraphics(f,sprintf('Fig3I_%01i.png',(i-1)*3+3),'Resolution',300,'BackgroundColor',[1 1 1]);
end

%% Fig 3 J

regions = [2,5,12];

allenHbTLow = mean(subAvg.Fig3.lowNE_FC_HbT,3);
allenHbTLow = allenHbTLow(regions,:);
allenHbTHigh = mean(subAvg.Fig3.highNE_FC_HbT,3);
allenHbTHigh = allenHbTHigh(regions,:);

for i = 1:3
    f = figure;f_plotAllenMap(allenHbTLow(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,sprintf('Fig3J_%01i.png',(i-1)*3+1),'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenHbTHigh(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,sprintf('Fig3J_%01i.png',(i-1)*3+2),'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenHbTHigh(i,:)-allenHbTLow(i,:),cmp=cmpbbr,mask=plotBM,cRange=0.1*[-1,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=0.7*[1 1 1]);
    colorbar off;exportgraphics(f,sprintf('Fig3J_%01i.png',(i-1)*3+3),'Resolution',300,'BackgroundColor',[1 1 1]);
end
