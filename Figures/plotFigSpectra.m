% calculate subject averages

NE_order = order(NE_Idx);

subAvg.Spectra.low_NE_Ca = NaN(300,numel(NE_order));
subAvg.Spectra.high_NE_Ca = NaN(300,numel(NE_order));
subAvg.Spectra.low_NE_HbT = NaN(300,numel(NE_order));
subAvg.Spectra.high_NE_HbT = NaN(300,numel(NE_order));
subAvg.Spectra.lowF_lowNE_Ca = NaN(12,numel(NE_order));
subAvg.Spectra.lowF_highNE_Ca = NaN(12,numel(NE_order));
subAvg.Spectra.medF_lowNE_Ca = NaN(12,numel(NE_order));
subAvg.Spectra.medF_highNE_Ca = NaN(12,numel(NE_order));
subAvg.Spectra.highF_lowNE_Ca = NaN(12,numel(NE_order));
subAvg.Spectra.highF_highNE_Ca = NaN(12,numel(NE_order));
subAvg.Spectra.lowF_lowNE_HbT = NaN(12,numel(NE_order));
subAvg.Spectra.lowF_highNE_HbT = NaN(12,numel(NE_order));
subAvg.Spectra.medF_lowNE_HbT = NaN(12,numel(NE_order));
subAvg.Spectra.medF_highNE_HbT = NaN(12,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.Spectra.low_NE_Ca(:,i) = mean(cat(2,spectra.low_NE_Ca{NE_order(i).Runs}),2);
    subAvg.Spectra.high_NE_Ca(:,i) = mean(cat(2,spectra.high_NE_Ca{NE_order(i).Runs}),2);
    subAvg.Spectra.low_NE_HbT(:,i) = mean(cat(2,spectra.low_NE_HbT{NE_order(i).Runs}),2);
    subAvg.Spectra.high_NE_HbT(:,i) = mean(cat(2,spectra.high_NE_HbT{NE_order(i).Runs}),2);
    subAvg.Spectra.lowF_lowNE_Ca(:,i) = mean(cat(2,spectra.lowF_lowNE_Ca{NE_order(i).Runs}),2);
    subAvg.Spectra.lowF_highNE_Ca(:,i) = mean(cat(2,spectra.lowF_highNE_Ca{NE_order(i).Runs}),2);
    subAvg.Spectra.medF_lowNE_Ca(:,i) = mean(cat(2,spectra.medF_lowNE_Ca{NE_order(i).Runs}),2);
    subAvg.Spectra.medF_highNE_Ca(:,i) = mean(cat(2,spectra.medF_highNE_Ca{NE_order(i).Runs}),2);
    subAvg.Spectra.highF_lowNE_Ca(:,i) = mean(cat(2,spectra.highF_lowNE_Ca{NE_order(i).Runs}),2);
    subAvg.Spectra.highF_highNE_Ca(:,i) = mean(cat(2,spectra.highF_highNE_Ca{NE_order(i).Runs}),2);
    subAvg.Spectra.lowF_lowNE_HbT(:,i) = mean(cat(2,spectra.lowF_lowNE_HbT{NE_order(i).Runs}),2);
    subAvg.Spectra.lowF_highNE_HbT(:,i) = mean(cat(2,spectra.lowF_highNE_HbT{NE_order(i).Runs}),2);
    subAvg.Spectra.medF_lowNE_HbT(:,i) = mean(cat(2,spectra.medF_lowNE_HbT{NE_order(i).Runs}),2);
    subAvg.Spectra.medF_highNE_HbT(:,i) = mean(cat(2,spectra.medF_highNE_HbT{NE_order(i).Runs}),2);
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

%% extract NE peaks

tBef = 50;
tAft = 50;

[peaks,trials] = f_findPeaks(spectra.NE,[0,0.02],10,0.005,[tBef,tAft]);

spg_trials_Ca = cell(numel(spectra.NE),1);
spg_trials_HbT = cell(numel(spectra.NE),1);
spg_trials_NE = cell(numel(spectra.NE),1);
spg_trials_pupil = cell(numel(spectra.NE),1);
subAvg.Spectra.spg_trials_Ca = NaN(tBef*10+tAft*10+1,300,12,numel(NE_order));
subAvg.Spectra.spg_trials_HbT = NaN(tBef*10+tAft*10+1,300,12,numel(NE_order));
subAvg.Spectra.spg_trials_NE = NaN(tBef*10+tAft*10+1,numel(NE_order));
subAvg.Spectra.spg_trials_pupil = NaN(tBef*10+tAft*10+1,numel(NE_order));

for i = 1:numel(spectra.NE)
    for idx = 1:numel(peaks{i})
        spg_trials_Ca{i}(:,:,:,idx) = spectra.SPG_Ca{i}(peaks{i}(idx)-tBef*10:peaks{i}(idx)+tAft*10,:,:);
        spg_trials_HbT{i}(:,:,:,idx) = spectra.SPG_HbT{i}(peaks{i}(idx)-tBef*10:peaks{i}(idx)+tAft*10,:,:);
        spg_trials_NE{i}(:,idx) = spectra.NE{i}(peaks{i}(idx)-tBef*10:peaks{i}(idx)+tAft*10);
        spg_trials_pupil{i}(:,idx) = Behavior.signals{i}(peaks{i}(idx)-tBef*10:peaks{i}(idx)+tAft*10,4);
    end
end

for i = 1:numel(NE_order)
    subAvg.Spectra.spg_trials_Ca(:,:,:,i) = mean(cat(4,spg_trials_Ca{NE_order(i).Runs}),4);
    subAvg.Spectra.spg_trials_HbT(:,:,:,i) = mean(cat(4,spg_trials_HbT{NE_order(i).Runs}),4);
    subAvg.Spectra.spg_trials_NE(:,i) = mean(cat(2,spg_trials_NE{NE_order(i).Runs}),2);
    subAvg.Spectra.spg_trials_pupil(:,i) = mean(cat(2,spg_trials_pupil{NE_order(i).Runs}),2);
end

%% Fig Spectra A
f = figure(Position=[100,100,1000,200]);
meanSig = mean(subAvg.Spectra.spg_trials_NE,2);
error = std(subAvg.Spectra.spg_trials_NE,0,2)/sqrt(M_NE);
f_plotLineError(-50:0.1:50,meanSig,error,color=c_GRAB,log=0,lineWidth=3);
axis off;

%% Fig Spectra B
f = figure(Position=[100,100,1000,300]);
imagesc(flipud(mean(subAvg.Spectra.spg_trials_Ca,[3,4])'));
axis off;
colormap cmpinf;
clim([0 0.15]);

%% Fig Spectra C
f = figure(Position=[100,100,1000,300]);
imagesc(flipud(mean(subAvg.Spectra.spg_trials_HbT,[3,4])'));
axis off;
colormap cmpinf;
clim([0 0.15]);

%% Fig Spectra D

f = figure;
meanSig = mean(subAvg.Spectra.low_NE_Ca,2);
error = std(subAvg.Spectra.low_NE_Ca,0,2)/sqrt(M_NE);
f_plotLineError(spectra.fr,spectra.fr.*meanSig,spectra.fr.*error,color=c_darkCyan,log=1,lineWidth=3);
meanSig = mean(subAvg.Spectra.high_NE_Ca,2);
error = std(subAvg.Spectra.high_NE_Ca,0,2)/sqrt(M_NE);
f_plotLineError(spectra.fr,spectra.fr.*meanSig,spectra.fr.*error,color=c_Orange,log=1,lineWidth=3);
set(gca,'YScale','linear','FontSize',14);
xlim([0.01 5]);
xlabel('F (Hz)');
ylabel('Power (normalized by F)');
legend('','low NE','','high NE');

f = figure;
meanSig = mean(subAvg.Spectra.low_NE_HbT,2);
error = std(subAvg.Spectra.low_NE_HbT,0,2)/sqrt(M_NE);
f_plotLineError(spectra.fr,spectra.fr.*meanSig,spectra.fr.*error,color=c_darkCyan,log=1,lineWidth=3);
meanSig = mean(subAvg.Spectra.high_NE_HbT,2);
error = std(subAvg.Spectra.high_NE_HbT,0,2)/sqrt(M_NE);
f_plotLineError(spectra.fr,spectra.fr.*meanSig,spectra.fr.*error,color=c_Orange,log=1,lineWidth=3);
set(gca,'YScale','linear','FontSize',14);
xlim([0.01 5]);
xlabel('F (Hz)');
ylabel('Power (normalized by F)');
legend('','low NE','','high NE');

%% Fig Spectra E

f = figure;f_plotAllenMap(mean(subAvg.Spectra.lowF_lowNE_Ca,2),cmp=cmpinf,mask=plotBM,cRange=[0,0.1]);
colorbar off;exportgraphics(f,'FigSpectraE_1.png','Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.Spectra.lowF_highNE_Ca,2),cmp=cmpinf,mask=plotBM,cRange=[0,0.1]);
colorbar off;exportgraphics(f,'FigSpectraE_2.png','Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.Spectra.lowF_highNE_Ca,2)-mean(subAvg.Spectra.lowF_lowNE_Ca,2),cmp=cmpbbr,mask=plotBM,cRange=0.03*[-1,1]);
colorbar off;exportgraphics(f,'FigSpectraE_3.png','Resolution',300,'BackgroundColor',[1 1 1]);
%%

f = figure;f_plotAllenMap(mean(subAvg.Spectra.medF_lowNE_Ca,2),cmp=cmpinf,mask=plotBM,cRange=[0,0.1]);
colorbar off;exportgraphics(f,'FigSpectraE_4.png','Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.Spectra.medF_highNE_Ca,2),cmp=cmpinf,mask=plotBM,cRange=[0,0.1]);
colorbar off;exportgraphics(f,'FigSpectraE_5.png','Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.Spectra.medF_highNE_Ca,2)-mean(subAvg.Spectra.medF_lowNE_Ca,2),cmp=cmpbbr,mask=plotBM,cRange=0.03*[-1,1]);
colorbar off;exportgraphics(f,'FigSpectraE_6.png','Resolution',300,'BackgroundColor',[1 1 1]);
%%

f = figure;f_plotAllenMap(mean(subAvg.Spectra.highF_lowNE_Ca,2),cmp=cmpinf,mask=plotBM,cRange=[0,0.05]);
colorbar off;exportgraphics(f,'FigSpectraE_7.png','Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.Spectra.highF_highNE_Ca,2),cmp=cmpinf,mask=plotBM,cRange=[0,0.05]);
colorbar off;exportgraphics(f,'FigSpectraE_8.png','Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.Spectra.highF_highNE_Ca,2)-mean(subAvg.Spectra.highF_lowNE_Ca,2),cmp=cmpbbr,mask=plotBM,cRange=0.02*[-1,1]);
colorbar off;exportgraphics(f,'FigSpectraE_9.png','Resolution',300,'BackgroundColor',[1 1 1]);
%% Fig Spectra F

f = figure;f_plotAllenMap(mean(subAvg.Spectra.lowF_lowNE_HbT,2),cmp=cmpinf,mask=plotBM,cRange=[0,0.05]);
colorbar off;exportgraphics(f,'FigSpectraF_1.png','Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.Spectra.lowF_highNE_HbT,2),cmp=cmpinf,mask=plotBM,cRange=[0,0.05]);
colorbar off;exportgraphics(f,'FigSpectraF_2.png','Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.Spectra.lowF_highNE_HbT,2)-mean(subAvg.Spectra.lowF_lowNE_HbT,2),cmp=cmpbbr,mask=plotBM,cRange=0.01*[-1,1]);
colorbar off;exportgraphics(f,'FigSpectraF_3.png','Resolution',300,'BackgroundColor',[1 1 1]);
%%

f = figure;f_plotAllenMap(mean(subAvg.Spectra.medF_lowNE_HbT,2),cmp=cmpinf,mask=plotBM,cRange=[0,0.1]);
colorbar off;exportgraphics(f,'FigSpectraF_4.png','Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.Spectra.medF_highNE_HbT,2),cmp=cmpinf,mask=plotBM,cRange=[0,0.1]);
colorbar off;exportgraphics(f,'FigSpectraF_5.png','Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.Spectra.medF_highNE_HbT,2)-mean(subAvg.Spectra.medF_lowNE_HbT,2),cmp=cmpbbr,mask=plotBM,cRange=0.02*[-1,1]);
colorbar off;exportgraphics(f,'FigSpectraF_6.png','Resolution',300,'BackgroundColor',[1 1 1]);

%%
function [peaks,trials] = f_findPeaks(NE,fwin,fr,th,t)

N = numel(NE);
peaks = cell(N,1);
trials = zeros(N,1);

for i = 1:N
    
    if isempty(NE{i})
        continue
    end

    sig = f_bpf(NE{i},fwin,fr);
    diffNE = diff(sig);
    transitions = diffNE>th;

    if sum(transitions)
        tIdx = find(diff(transitions)==1);
        if isempty(tIdx)
            continue
        end
        tDown = find(diff(transitions)==-1);
        tDown(tDown<tIdx(1)) = [];
        if numel(tIdx)>numel(tDown)
            tIdx(numel(tDown)+1:end) = [];
        end

        transitions = round(mean([tIdx tDown],2));
        transitions(transitions<t(1)*10+1) = [];
        transitions(transitions>numel(sig)-t(2)*10) = [];
        peaks{i} = transitions;
        trials(i) = numel(transitions);
    end
end

end

