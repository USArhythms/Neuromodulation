% calculate subject averages

NE_order = order(NE_Idx);

subAvg.FC.lowF_lowNE_Ca = NaN(12,12,numel(NE_order));
subAvg.FC.lowF_highNE_Ca = NaN(12,12,numel(NE_order));
subAvg.FC.medF_lowNE_Ca = NaN(12,12,numel(NE_order));
subAvg.FC.medF_highNE_Ca = NaN(12,12,numel(NE_order));
subAvg.FC.highF_lowNE_Ca = NaN(12,12,numel(NE_order));
subAvg.FC.highF_highNE_Ca = NaN(12,12,numel(NE_order));
subAvg.FC.lowF_lowNE_HbT = NaN(12,12,numel(NE_order));
subAvg.FC.lowF_highNE_HbT = NaN(12,12,numel(NE_order));
subAvg.FC.medF_lowNE_HbT = NaN(12,12,numel(NE_order));
subAvg.FC.medF_highNE_HbT = NaN(12,12,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.FC.lowF_lowNE_Ca(:,:,i) = mean(cat(3,FC_fr.lowF_lowNE_Ca{NE_order(i).Runs}),3);
    subAvg.FC.lowF_highNE_Ca(:,:,i) = mean(cat(3,FC_fr.lowF_highNE_Ca{NE_order(i).Runs}),3);
    subAvg.FC.medF_lowNE_Ca(:,:,i) = mean(cat(3,FC_fr.medF_lowNE_Ca{NE_order(i).Runs}),3);
    subAvg.FC.medF_highNE_Ca(:,:,i) = mean(cat(3,FC_fr.medF_highNE_Ca{NE_order(i).Runs}),3);
    subAvg.FC.highF_lowNE_Ca(:,:,i) = mean(cat(3,FC_fr.highF_lowNE_Ca{NE_order(i).Runs}),3);
    subAvg.FC.highF_highNE_Ca(:,:,i) = mean(cat(3,FC_fr.highF_highNE_Ca{NE_order(i).Runs}),3);
    subAvg.FC.lowF_lowNE_HbT(:,:,i) = mean(cat(3,FC_fr.lowF_lowNE_HbT{NE_order(i).Runs}),3);
    subAvg.FC.lowF_highNE_HbT(:,:,i) = mean(cat(3,FC_fr.lowF_highNE_HbT{NE_order(i).Runs}),3);
    subAvg.FC.medF_lowNE_HbT(:,:,i) = mean(cat(3,FC_fr.medF_lowNE_HbT{NE_order(i).Runs}),3);
    subAvg.FC.medF_highNE_HbT(:,:,i) = mean(cat(3,FC_fr.medF_highNE_HbT{NE_order(i).Runs}),3);
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

%% Fig FC A

figName = 'FigFcA';

f = figure;
f_plotFC(mean(subAvg.FC.lowF_lowNE_Ca,3),1,cmp=cmpvir,bounds=[0 1],title='Low NE Ca++ Connectivity',clabel='r');
title('');colorbar off;
exportgraphics(f,[figName '_1.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;
f_plotFC(mean(subAvg.FC.lowF_highNE_Ca,3),1,cmp=cmpvir,bounds=[0 1],title='High NE Ca++ Connectivity',clabel='r');
title('');colorbar off;
exportgraphics(f,[figName '_2.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,FC_fr.lowF_lowNE_Ca{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,FC_fr.lowF_highNE_Ca{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

f = figure;
f_plotFC(mean(subAvg.FC.lowF_highNE_Ca,3)-mean(subAvg.FC.lowF_lowNE_Ca,3),0,cmp=cmpbbr,bounds=0.3*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));
title('');colorbar off;
exportgraphics(f,[figName '_3.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig FC B

figName = 'FigFcB';

f = figure;
f_plotFC(mean(subAvg.FC.medF_lowNE_Ca,3),1,cmp=cmpvir,bounds=[0 1],title='Low NE Ca++ Connectivity',clabel='r');
title('');colorbar off;
exportgraphics(f,[figName '_1.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;
f_plotFC(mean(subAvg.FC.medF_highNE_Ca,3),1,cmp=cmpvir,bounds=[0 1],title='High NE Ca++ Connectivity',clabel='r');
title('');colorbar off;
exportgraphics(f,[figName '_2.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,FC_fr.medF_lowNE_Ca{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,FC_fr.medF_highNE_Ca{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

f = figure;
f_plotFC(mean(subAvg.FC.medF_highNE_Ca,3)-mean(subAvg.FC.medF_lowNE_Ca,3),0,cmp=cmpbbr,bounds=0.2*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));
title('');colorbar off;
exportgraphics(f,[figName '_3.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig FC C

figName = 'FigFcC';

f = figure;
f_plotFC(mean(subAvg.FC.highF_lowNE_Ca,3),1,cmp=cmpvir,bounds=[0 1],title='Low NE Ca++ Connectivity',clabel='r');
title('');colorbar off;
exportgraphics(f,[figName '_1.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;
f_plotFC(mean(subAvg.FC.highF_highNE_Ca,3),1,cmp=cmpvir,bounds=[0 1],title='High NE Ca++ Connectivity',clabel='r');
title('');colorbar off;
exportgraphics(f,[figName '_2.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,FC_fr.highF_lowNE_Ca{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,FC_fr.highF_highNE_Ca{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

f = figure;
f_plotFC(mean(subAvg.FC.highF_highNE_Ca,3)-mean(subAvg.FC.highF_lowNE_Ca,3),0,cmp=cmpbbr,bounds=0.1*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));
title('');colorbar off;
exportgraphics(f,[figName '_3.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig FC D

figName = 'FigFcD';

f = figure;
f_plotFC(mean(subAvg.FC.lowF_lowNE_HbT,3),1,cmp=cmpvir,bounds=[0 1],title='Low NE Ca++ Connectivity',clabel='r');
title('');colorbar off;
exportgraphics(f,[figName '_1.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;
f_plotFC(mean(subAvg.FC.lowF_highNE_HbT,3),1,cmp=cmpvir,bounds=[0 1],title='High NE Ca++ Connectivity',clabel='r');
title('');colorbar off;
exportgraphics(f,[figName '_2.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,FC_fr.lowF_lowNE_HbT{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,FC_fr.lowF_highNE_HbT{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

f = figure;
f_plotFC(mean(subAvg.FC.lowF_highNE_HbT,3)-mean(subAvg.FC.lowF_lowNE_HbT,3),0,cmp=cmpbbr,bounds=0.15*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));
title('');colorbar off;
exportgraphics(f,[figName '_3.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig FC E

figName = 'FigFcE';

f = figure;
f_plotFC(mean(subAvg.FC.medF_lowNE_HbT,3),1,cmp=cmpvir,bounds=[0 1],title='Low NE Ca++ Connectivity',clabel='r');
title('');colorbar off;
exportgraphics(f,[figName '_1.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;
f_plotFC(mean(subAvg.FC.medF_highNE_HbT,3),1,cmp=cmpvir,bounds=[0 1],title='High NE Ca++ Connectivity',clabel='r');
title('');colorbar off;
exportgraphics(f,[figName '_2.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,FC_fr.medF_lowNE_HbT{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,FC_fr.medF_highNE_HbT{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

f = figure;
f_plotFC(mean(subAvg.FC.medF_highNE_HbT,3)-mean(subAvg.FC.medF_lowNE_HbT,3),0,cmp=cmpbbr,bounds=0.1*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));
title('');colorbar off;
exportgraphics(f,[figName '_3.png'],'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig FC F

regions = [2,5,12];

allenCaLow = mean(subAvg.FC.lowF_lowNE_Ca,3);
allenCaLow = allenCaLow(regions,:);
allenCaHigh = mean(subAvg.FC.lowF_highNE_Ca,3);
allenCaHigh = allenCaHigh(regions,:);

figName = 'FigFcF';

for i = 1:3
    f = figure;f_plotAllenMap(allenCaLow(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+1)],'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenCaHigh(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+2)],'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenCaHigh(i,:)-allenCaLow(i,:),cmp=cmpbbr,mask=plotBM,cRange=0.3*[-1,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=0.7*[1 1 1]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+3)],'Resolution',300,'BackgroundColor',[1 1 1]);
end

%% Fig FC G

regions = [2,5,12];

allenCaLow = mean(subAvg.FC.medF_lowNE_Ca,3);
allenCaLow = allenCaLow(regions,:);
allenCaHigh = mean(subAvg.FC.medF_highNE_Ca,3);
allenCaHigh = allenCaHigh(regions,:);

figName = 'FigFcG';

for i = 1:3
    f = figure;f_plotAllenMap(allenCaLow(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+1)],'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenCaHigh(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+2)],'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenCaHigh(i,:)-allenCaLow(i,:),cmp=cmpbbr,mask=plotBM,cRange=0.2*[-1,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=0.7*[1 1 1]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+3)],'Resolution',300,'BackgroundColor',[1 1 1]);
end

%% Fig FC H

regions = [2,5,12];

allenCaLow = mean(subAvg.FC.highF_lowNE_Ca,3);
allenCaLow = allenCaLow(regions,:);
allenCaHigh = mean(subAvg.FC.highF_highNE_Ca,3);
allenCaHigh = allenCaHigh(regions,:);

figName = 'FigFcH';

for i = 1:3
    f = figure;f_plotAllenMap(allenCaLow(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+1)],'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenCaHigh(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+2)],'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenCaHigh(i,:)-allenCaLow(i,:),cmp=cmpbbr,mask=plotBM,cRange=0.1*[-1,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=0.7*[1 1 1]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+3)],'Resolution',300,'BackgroundColor',[1 1 1]);
end

%% Fig FC I

regions = [2,5,12];

allenCaLow = mean(subAvg.FC.lowF_lowNE_HbT,3);
allenCaLow = allenCaLow(regions,:);
allenCaHigh = mean(subAvg.FC.lowF_highNE_HbT,3);
allenCaHigh = allenCaHigh(regions,:);

figName = 'FigFcI';

for i = 1:3
    f = figure;f_plotAllenMap(allenCaLow(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+1)],'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenCaHigh(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+2)],'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenCaHigh(i,:)-allenCaLow(i,:),cmp=cmpbbr,mask=plotBM,cRange=0.15*[-1,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=0.7*[1 1 1]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+3)],'Resolution',300,'BackgroundColor',[1 1 1]);
end

%% Fig FC J

regions = [2,5,12];

allenCaLow = mean(subAvg.FC.medF_lowNE_HbT,3);
allenCaLow = allenCaLow(regions,:);
allenCaHigh = mean(subAvg.FC.medF_highNE_HbT,3);
allenCaHigh = allenCaHigh(regions,:);

figName = 'FigFcJ';

for i = 1:3
    f = figure;f_plotAllenMap(allenCaLow(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+1)],'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenCaHigh(i,:),cmp=cmpvir,mask=plotBM,cRange=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+2)],'Resolution',300,'BackgroundColor',[1 1 1]);
    f = figure;f_plotAllenMap(allenCaHigh(i,:)-allenCaLow(i,:),cmp=cmpbbr,mask=plotBM,cRange=0.1*[-1,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=0.7*[1 1 1]);
    colorbar off;exportgraphics(f,[figName sprintf('_%01i.png',(i-1)*3+3)],'Resolution',300,'BackgroundColor',[1 1 1]);
end

