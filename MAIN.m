%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          Main Processing Script for Neuromodulation Analysis
%
% Reads nwb files located at https://dandiarchive.org/dandiset/001211

%% add path to local matNWB installation and child directories

% ignore if matNWB is installed via MATLAB
addpath('/projectnb/devorlab/bcraus/AnalysisCode/new_processing/matnwb');

addPaths(fileparts(mfilename('fullpath')));

%% load nwb file and extract relevant variables

% local path to nwb file directory (load from DANDI archives)
nwbDir = '/projectnb/devorlab/bcraus/AnalysisCode/NWB/nwb_files_HD';
nwb_list = f_sortNWB(nwbDir);

% temporary index to analyze only a single run
nwbI = 90;

nwb = nwbRead(nwb_list(nwbI).Path);

[rfp,gfp,rfp_HD,gfp_HD,Hb,HbO,HbT,Whisking,Pupil,Accelerometer,brain_mask,vessel_mask,allen_masks,fs,mouseInfo,sessionInfo,behCam] = f_extractNWB(nwb);

GRAB = mouseInfo.GRAB;

%% perform analysis for Fig 1

Fig1 = struct;

ds = 32; % downsampling factor for estimation of initial timing parameters of IRF
numThreads = 4; % number of cores to use in optimization algorithm
win = [0 10]; % IRF kernel range (s)
corrWin = [15 3]; % window for temporal IRF performance analysis (s)

Fig1.irf_win = win;
Fig1.perf_dt_win = corrWin;

% low pass filter HbT below 0.5 Hz, for unfiltered use HbT
HbT_low = f_bpf(HbT,[0, 0.5],fs,3);

% estimate invariant single IRF models

Fig1.IRFx1_inv = struct;

[Fig1.IRFx1_inv.perf,Fig1.IRFx1_inv.IRF,Fig1.IRFx1_inv.params,tmpCorr] = f_1xIRF(HbT_low,rfp_HD,win,fs,brain_mask,4,corrWin*fs,numThreads);
Fig1.IRFx1_inv.perf_dt = squeeze(mean(tmpCorr.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));

% estimate medial SSp single IRF model 

Fig1.IRFx1_SSp = struct;

SSp = sum(allen_masks(:,:,[4 5]),3,'omitnan');
SSp(SSp==0) = NaN;

[Fig1.IRFx1_SSp.perf,Fig1.IRFx1_SSp.IRF,Fig1.IRFx1_SSp.params,tmpCorr] = f_1xIRF(HbT_low,rfp_HD,win,fs,SSp.*brain_mask,4,corrWin*fs,numThreads);
Fig1.IRFx1_SSp.perf_dt = squeeze(mean(tmpCorr.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));

% estimate variant IRF model

Fig1.IRFx1_var = struct;

[Fig1.IRFx1_var.perf,Fig1.IRFx1_var.IRF,Fig1.IRFx1_var.params,tmpCorr] = f_1xIRF_varWeights(HbT_low,rfp_HD,win,fs,brain_mask,ds,corrWin*fs,numThreads);
Fig1.IRFx1_var.perf_dt = squeeze(mean(tmpCorr.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
Fig1.IRFx1_var.IRF_allen = squeeze(mean(Fig1.IRFx1_var.IRF.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));

% apply moving window to GRAB and pupil

Fig1.GRAB = squeeze(mean(gfp_HD.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
Fig1.GRAB = movmean(Fig1.GRAB,corrWin(2)*fs);
Fig1.GRAB = Fig1.GRAB(corrWin(1)*fs/2:corrWin(2)*fs:size(Fig1.GRAB,1)-corrWin(1)*fs/2,:);

Fig1.Pupil = movmean(Pupil,corrWin(2)*fs);
Fig1.Pupil = Fig1.Pupil(corrWin(1)*fs/2:corrWin(2)*fs:size(Fig1.Pupil,1)-corrWin(1)*fs/2,:);

% calculate correlation coefficient - GRAB and pupil vs. SSp IRF performance
Fig1.SSp_perf_vs_pupil = corrcoef([Fig1.Pupil, Fig1.IRFx1_SSp.perf_dt]);
Fig1.SSp_perf_vs_pupil = Fig1.SSp_perf_vs_pupil(2:13)';

Fig1.SSp_perf_vs_GRAB = f_corr(Fig1.GRAB, Fig1.IRFx1_SSp.perf_dt, 1)';

%% perform analysis for Fig 2

if string(GRAB) == "GRAB_NE"
    Fig2 = struct;
    
    ds = 32; % downsampling factor for estimation of initial timing parameters of IRF
    numThreads = 4; % number of cores to use in optimization algorithm
    win = [-5 10]; % IRF kernel range (s)
    
    Fig2.irf_win = win;
    
    HbT_low = f_bpf(HbT,[0, 0.5],fs,3);
    
    % estimate linear regression model
    
    Fig2.LR = struct;
    [Fig2.LR.perf,Fig2.LR.params] = f_LR_varWeights(HbT_low,rfp_HD,gfp_HD,win,fs,brain_mask,ds,numThreads);
    
    % estimate double IRF model - varying weights
    
    Fig2.IRFx2 = struct;
    [Fig2.IRFx2.perf,Fig2.IRFx2.IRF,Fig2.IRFx2.params] = f_2alphaDeconvolve(HbT_low,rfp_HD,gfp_HD,win,fs,brain_mask,ds,numThreads);
    Fig2.IRFx2.params.A = Fig2.IRFx2.params.A*sum(Fig2.IRFx2.IRF(:,1));
    Fig2.IRFx2.params.B = Fig2.IRFx2.params.B*sum(Fig2.IRFx2.IRF(:,2));

end

%% perform analysis for Fig 3

Fig3 = struct;

win = [30 6]; % length of sliding connectivity window (s)
Fig3.win = win;

% extract Ca and HbT ROIs

HbT_low = f_bpf(HbT,[0, 0.5],fs,3);

Fig3.Ca = squeeze(mean(rfp_HD./std(rfp_HD,0,3).*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
Fig3.HbT = squeeze(mean(HbT_low./std(HbT_low,0,3).*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));

% extract hemisphere-wide GRAB signal
Fig3.GRAB = squeeze(mean(gfp_HD.*brain_mask,[1, 2],'omitnan'));
Fig3.GRAB = movmean(Fig3.GRAB,win(2)*fs);
Fig3.GRAB = Fig3.GRAB(win(1)*fs/2:win(2)*fs:size(Fig3.GRAB,1)-win(1)*fs/2);

% percentile GRAB signal
bounds = prctile(Fig3.GRAB,[30, 70]);

lowNE = Fig3.GRAB < bounds(1);
highNE = Fig3.GRAB > bounds(2);

% calculate dynamic FC for Ca and HbT
Fig3.FC.Ca = f_funConGram(Fig3.Ca,win*fs);
Fig3.FC.HbT = f_funConGram(Fig3.HbT,win*fs);

% average dynamic FC for low and high NE
Fig3.lowNE_FC.Ca = mean(Fig3.FC.Ca(:,:,lowNE),3);
Fig3.lowNE_FC.HbT = mean(Fig3.FC.HbT(:,:,lowNE),3);

Fig3.highNE_FC.Ca = mean(Fig3.FC.Ca(:,:,highNE),3);
Fig3.highNE_FC.HbT = mean(Fig3.FC.HbT(:,:,highNE),3);

% calculate correlation between Ca and HbT FC
trilIdx = ~tril(ones(12));

N = size(Fig3.FC.Ca,3);

tmpCa = reshape(Fig3.FC.Ca,[],N);
tmpHbT = reshape(Fig3.FC.HbT,[],N);

Fig3.FC.Ca_vs_HbT = f_corr(tmpCa(trilIdx,:),tmpHbT(trilIdx,:),1)';

% calculate correlation measures - FC Ca, FC HbT, and FC Ca/HbT vs. GRAB

Fig3.r.FC_Ca_vs_GRAB = f_corr(Fig3.FC.Ca,permute(Fig3.GRAB,[2, 3, 1]),3);
Fig3.r.FC_HbT_vs_GRAB = f_corr(Fig3.FC.HbT,permute(Fig3.GRAB,[2, 3, 1]),3);
Fig3.r.FC_Ca_HbT_vs_GRAB = f_corr(Fig3.FC.Ca_vs_HbT,Fig3.GRAB,1);

%% analysis for behavior supplementary

[spectra,fr] = f_hemSpectra(rfp_HD,10,[0 5],brain_mask,[5,9]);
[C,phi,fr] = f_hemCoherence(rfp_HD,HbT,10,brain_mask,[5,9]);

%% organize main figures

NVC_Figures = struct;
NVC_Figures.Mouse = mouseInfo.ID;
NVC_Figures.Date = sessionInfo.Date;
NVC_Figures.Run = sessionInfo.Run;
NVC_Figures.iRun = sessionInfo.interRun;

% organize Figure 1
NVC_Figures.Fig1C.IRF = Fig1.IRFx1_inv.IRF;
NVC_Figures.Fig1C.performance = Fig1.IRFx1_inv.perf;

NVC_Figures.Fig1D.IRF = Fig1.IRFx1_SSp.IRF;
NVC_Figures.Fig1D.performance = Fig1.IRFx1_SSp.perf;

NVC_Figures.Fig1E.IRF_SSp = mean(Fig1.IRFx1_var.IRF_allen(:,4:5),2);
NVC_Figures.Fig1E.IRF_MOs = Fig1.IRFx1_var.IRF_allen(:,2);
NVC_Figures.Fig1E.performance = Fig1.IRFx1_var.perf;

NVC_Figures.Fig1F.NE = Fig1.GRAB(:,2);
NVC_Figures.Fig1F.performance_dt = Fig1.IRFx1_SSp.perf_dt(:,2);
NVC_Figures.Fig1F.R = Fig1.SSp_perf_vs_GRAB(2);

NVC_Figures.Fig1G.R = Fig1.SSp_perf_vs_GRAB;

% organize Figure 2
NVC_Figures.Fig2B.A = Fig2.LR.params.A;
NVC_Figures.Fig2B.B = Fig2.LR.params.B;
NVC_Figures.Fig2C.tA = Fig2.LR.params.tA;
NVC_Figures.Fig2C.tB = Fig2.LR.params.tB;
NVC_Figures.Fig2D.performance = Fig2.LR.perf;
NVC_Figures.Fig2E.diff_performance = Fig2.LR.perf-Fig1.IRFx1_inv.perf;
NVC_Figures.Fig2F.A = Fig2.IRFx2.params.A;
NVC_Figures.Fig2F.B = Fig2.IRFx2.params.B;
NVC_Figures.Fig2G.IRF_Ca = Fig2.IRFx2.IRF(:,1);
NVC_Figures.Fig2G.IRF_NE = Fig2.IRFx2.IRF(:,2);
NVC_Figures.Fig2H.performance = Fig2.IRFx2.perf;
NVC_Figures.Fig2I.diff_performance = Fig2.IRFx2.perf-Fig1.IRFx1_inv.perf;

% organize Figure 3
NVC_Figures.Fig3D.NE = Fig3.GRAB;
NVC_Figures.Fig3D.FC_Ca = squeeze(Fig3.FC.Ca(2,5,:));
NVC_Figures.Fig3D.FC_HbT = squeeze(Fig3.FC.HbT(2,5,:));
NVC_Figures.Fig3E.NE = Fig3.GRAB;
NVC_Figures.Fig3E.R = Fig3.FC.Ca_vs_HbT;
NVC_Figures.Fig3F.R = [Fig3.r.FC_Ca_vs_GRAB(2,5), Fig3.r.FC_HbT_vs_GRAB(2,5), Fig3.r.FC_Ca_HbT_vs_GRAB];
NVC_Figures.Fig3G.low = Fig3.lowNE_FC.Ca;
NVC_Figures.Fig3G.high = Fig3.highNE_FC.Ca;
NVC_Figures.Fig3G.diff = Fig3.highNE_FC.Ca-Fig3.lowNE_FC.Ca;
NVC_Figures.Fig3H.low = Fig3.lowNE_FC.HbT;
NVC_Figures.Fig3H.high = Fig3.highNE_FC.HbT;
NVC_Figures.Fig3H.diff = Fig3.highNE_FC.HbT-Fig3.lowNE_FC.HbT;

%% additional functions

function addPaths(mainDir)
    addDirs = genpath(mainDir);
    addDirs = strsplit(addDirs,pathsep);
    addDirs = addDirs(~cellfun(@(x) any(startsWith(strsplit(x,filesep), '.')), addDirs));
    addDirs = strjoin(addDirs,pathsep);
    addpath(addDirs);
end
