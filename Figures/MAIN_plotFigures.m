%% organize data

[log, order, settings, Fig1, Fig2, Fig3, Behavior, Hb_model, unfiltered, shuffled, NE_reg, spectra, FC_fr] = f_organizeData('/projectnb/devorlab/bcraus/AnalysisCode/Neuromodulation/Analysis');
load('/projectnb/devorlab/bcraus/refAllen.mat');

M = numel(order);

BM = f_regImages(settings.brain_mask,refParcellation,settings.parcellation,1);
NE_Idx = strcmp({order.GRAB},'GRAB_NE');
ACh_Idx = strcmp({order.GRAB},'GRAB_ACh');

%% plot Main figures
plotFig1
plotFig2

%% plot Supplementary figures

