%% organize data

[log, order, settings, Fig1, Fig2, Fig3, Behavior, Hb_model, unfiltered, shuffled, NE_reg, spectra, FC_fr, GRAB_FC] = f_organizeData('/projectnb/devorlab/bcraus/AnalysisCode/Neuromodulation/Analysis');

parentDir = f_path();
load(fullfile(parentDir,'Figures/plot_types/refAllen.mat'));

M = numel(order);

BM = f_regImages(settings.brain_mask,refParcellation,settings.parcellation,1);
NE_Idx = strcmp({order.GRAB},'GRAB_NE');
ACh_Idx = strcmp({order.GRAB},'GRAB_ACh');

set(0,'defaultfigurecolor',[1 1 1]);

M_NE = sum(NE_Idx);
M_ACh = sum(ACh_Idx);

%% plot Main figures
plotFig1
plotFig2

%% plot Supplementary figures

