%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%          Main Processing Script for Neuromodulation Analysis
% 
% Reads nwb files located at https://dandiarchive.org/dandiset/001543 and
% performs IRF modeling, functional connectivity, and spectral analysis.

%% add path to local matNWB installation and child directories

% get path to Neuromodulation directory
parentDir = f_path();
addPaths(parentDir);

% add path to matNWB directory - ignore if matNWB is installed
addpath('/projectnb/devorlab/bcraus/AnalysisCode/new_processing/matnwb');

%% load nwb file and extract relevant variables

% local path to nwb file directory (load from DANDI archives)
nwbDir = fullfile(parentDir,'data');
nwb_list = f_sortNWB(nwbDir);

%% define files and input parameters

saveMetadata = false;

files = struct;
files.metadata = fullfile(parentDir,'Analysis');
files.save = fullfile(parentDir,'results');

%%
% loop through all .nwb files in nwb_list
for nwbI = 1:numel(nwb_list)

    fprintf(['Analyzing ', nwb_list(nwbI).Path, '...']);
    
    % read and load data from nwb file
    nwb = nwbRead(nwb_list(nwbI).Path);
    [~,gfp,rfp_HD,gfp_HD,Hb,HbO,HbT,Whisking,Pupil,Accelerometer,...
        brain_mask,vessel_mask,allen_masks,fs,mouseInfo,sessionInfo]...
        = f_extractNWB(nwb);
        
    %% perform analysis for Fig 1
    
    fprintf('\n\tFigure 1 Analysis...');
    
    Fig1 = struct;
    
    ds = 32; % downsampling factor for estimation of initial timing parameters of IRF
    numThreads = 8; % number of cores to use in optimization algorithm
    win = [0 10]; % IRF kernel range (s)
    corrWin = [15 3]; % window for temporal IRF performance analysis (s)
    
    Fig1.irf_win = win;
    Fig1.perf_dt_win = corrWin;
    
    % low pass filter HbT below 0.5 Hz, for unfiltered use HbT
    HbT_low = f_bpf(HbT,[0,0.5],fs,3);
    
    Fig1.Ca_allen = f_parcellate(rfp_HD,allen_masks);
    Fig1.GRAB_allen = f_parcellate(gfp_HD,allen_masks);
    Fig1.HbT_allen = f_parcellate(HbT,allen_masks);
    Fig1.Pupil_full = Pupil;
    Fig1.Whisking = Whisking;
    Fig1.Accelerometer = Accelerometer;
    Fig1.GRAB_type = mouseInfo.GRAB;
    Fig1.allen_masks = allen_masks;

    % estimate invariant single IRF models
    [Fig1.IRFx1_inv.perf,Fig1.IRFx1_inv.IRF,Fig1.IRFx1_inv.params,tmpCorr] = f_1xIRF(HbT_low,rfp_HD,win,fs,brain_mask,4,corrWin*fs,numThreads);
    Fig1.IRFx1_inv.perf_dt = squeeze(mean(tmpCorr.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
    
    % estimate medial SSp single IRF model 
    SSp = sum(allen_masks(:,:,[4 5]),3,'omitnan');
    SSp(SSp==0) = NaN;
    
    [Fig1.IRFx1_SSp.perf,Fig1.IRFx1_SSp.IRF,Fig1.IRFx1_SSp.params,tmpCorr] = f_1xIRF(HbT_low,rfp_HD,win,fs,SSp.*brain_mask,4,corrWin*fs,numThreads);
    Fig1.IRFx1_SSp.perf_dt = squeeze(mean(tmpCorr.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
    
    % estimate variant IRF model
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

    %% GRAB connectivity

    GRAB_FC = struct;

    GRAB_FC.GRAB = squeeze(mean(gfp_HD.*brain_mask.*vessel_mask.*permute(allen_masks,[1 2 4 3]),[1,2],'omitnan'));
    GRAB_FC.GRAB_global = squeeze(mean(gfp_HD.*brain_mask.*vessel_mask,[1,2],'omitnan'));
    GRAB_FC.GRAB_norm = squeeze(mean(gfp_HD./std(gfp_HD,0,3).*brain_mask.*vessel_mask,[1,2],'omitnan'));
    GRAB_FC.FC = corrcoef(GRAB_FC.GRAB);
    GRAB_FC.FC_detrend = corrcoef(detrend(GRAB_FC.GRAB));

    %% perform analysis for Fig 2
    
    if string(mouseInfo.GRAB) == "GRAB_NE"
        fprintf('\n\tFigure 2 Analysis...');
        Fig2 = struct;
        
        ds = 32; % downsampling factor for estimation of initial timing parameters of IRF
        numThreads = 4; % number of cores to use in optimization algorithm
        win = [-5 10]; % IRF kernel range (s)
        
        Fig2.irf_win = win;
        Fig2.allen_masks = allen_masks;

        HbT_low = f_bpf(HbT,[0, 0.5],fs,3);
        
        globalNE = f_parcellate(gfp_HD./std(gfp_HD,0,3),brain_mask);
        globalNE = permute(globalNE,[2,3,1]);
        Fig2.Ca_norm = f_parcellate(rfp_HD./std(rfp_HD,0,3),allen_masks);
        Fig2.HbT_norm = f_parcellate(HbT_low./std(HbT_low,0,3),allen_masks);

        Fig2.LR = struct;
        Fig2.NE = squeeze(globalNE);

        % estimate linear regression model

        [Fig2.LR.perf,Fig2.LR.params,Fig2.LR_Ca,Fig2.LR_NE] = f_LR_varWeights(HbT_low,rfp_HD,globalNE.*ones(size(rfp_HD,[1,2])),win,fs,brain_mask,ds,numThreads);
        Fig2.LR_Ca = f_parcellate(Fig2.LR_Ca,brain_mask.*allen_masks);
        Fig2.LR_NE = f_parcellate(Fig2.LR_NE,brain_mask.*allen_masks);

        % estimate double IRF model - varying weights
        
        Fig2.IRFx2 = struct;
        [Fig2.IRFx2.perf,Fig2.IRFx2.IRF,Fig2.IRFx2.params,Fig2.IRFx2_Ca,Fig2.IRFx2_NE] = f_2alphaDeconvolve(HbT_low,rfp_HD,globalNE.*ones(size(rfp_HD,[1,2])),win,fs,brain_mask,ds,numThreads);
        Fig2.IRFx2.params.A = Fig2.IRFx2.params.A*sum(Fig2.IRFx2.IRF(:,1));
        Fig2.IRFx2.params.B = Fig2.IRFx2.params.B*sum(Fig2.IRFx2.IRF(:,2));
        Fig2.IRFx2_Ca = f_parcellate(Fig2.IRFx2_Ca,brain_mask.*allen_masks);
        Fig2.IRFx2_NE = f_parcellate(Fig2.IRFx2_NE,brain_mask.*allen_masks);

    end

    %% perform analysis for Fig 3
    
    if string(mouseInfo.GRAB) == "GRAB_NE"
        fprintf('\n\tFigure 3 Analysis...');
        Fig3 = struct;
        
        win = [30 6]; % length of sliding connectivity window (s)
        Fig3.win = win;
        
        % extract Ca and HbT ROIs
        
        HbT_low = f_bpf(HbT,[0, 0.5],fs,3);
        
        Fig3.Pupil = Pupil;
        Fig3.Whisking = Whisking;
        Fig3.Accelerometer = Accelerometer;

        Fig3.Ca = squeeze(mean(rfp_HD./std(rfp_HD,0,3).*vessel_mask.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
        Fig3.HbT = squeeze(mean(HbT_low./std(HbT_low,0,3).*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
        
        % extract hemisphere-wide GRAB signal
        Fig3.GRAB = squeeze(mean(gfp_HD.*vessel_mask,[1, 2],'omitnan'));
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
    end

    %% plot Fig 3
    
    if string(mouseInfo.GRAB) == "GRAB_NE"
        plotFig3(Fig3,files.save);
    end
    
    %% analysis for behavior supplementary
    fprintf('\n\tBehavior Analysis...');

    Behavior = struct;
    Behavior.tapers = [5 9];
    
    comb_mask = brain_mask.*vessel_mask;
    
    % spectrum
    [Behavior.SPC.rfp_HD,Behavior.SPC.fr] = f_hemSpectra(rfp_HD,fs,[0 round(fs/2)],comb_mask,Behavior.tapers);
    [Behavior.SPC.gfp_HD,Behavior.SPC.fr] = f_hemSpectra(gfp_HD,fs,[0 round(fs/2)],comb_mask,Behavior.tapers);
    [Behavior.SPC.gfp,Behavior.SPC.fr] = f_hemSpectra(gfp,fs,[0 round(fs/2)],comb_mask,Behavior.tapers);
    [Behavior.SPC.HbT,Behavior.SPC.fr] = f_hemSpectra(HbT,fs,[0 round(fs/2)],brain_mask,Behavior.tapers);
    
    % coherence
    [Behavior.COH.rfp_HD_HbT,Behavior.PHI.rfp_HD_HbT] = f_hemCoherence(rfp_HD,HbT,fs,comb_mask,Behavior.tapers);
    [Behavior.COH.gfp_HD_HbT,Behavior.PHI.gfp_HD_HbT] = f_hemCoherence(gfp_HD,HbT,fs,comb_mask,Behavior.tapers);
    [Behavior.COH.gfp_HbT,Behavior.PHI.gfp_HbT] = f_hemCoherence(gfp,HbT,fs,comb_mask,Behavior.tapers);
    [Behavior.COH.rfp_HD_gfp_HD,Behavior.PHI.rfp_HD_gfp_HD] = f_hemCoherence(rfp_HD,gfp_HD,fs,comb_mask,Behavior.tapers);
    
    % xcorr
    [Behavior.XC.rfp_HD_HbT,Behavior.XC.lag] = f_hemLag_dT(rfp_HD,HbT_low,fs,[-10 10],comb_mask);
    [Behavior.XC.gfp_HD_HbT,Behavior.XC.lag] = f_hemLag_dT(gfp_HD,HbT_low,fs,[-10 10],comb_mask);
    [Behavior.XC.gfp_HbT,Behavior.XC.lag] = f_hemLag_dT(gfp,HbT_low,fs,[-10 10],comb_mask);
    [Behavior.XC.rfp_HD_gfp_HD,Behavior.XC.lag] = f_hemLag_dT(rfp_HD,gfp_HD,fs,[-10 10],comb_mask);
    
    % correlation
    Behavior.R.rfp_HD_low_gfp_HD_low = f_corr(f_bpf(rfp_HD,[0, 0.5],fs,3),f_bpf(gfp_HD,[0, 0.5],fs,3),3);
    Behavior.R.rfp_HD_low_HbT_low = f_corr(f_bpf(rfp_HD,[0, 0.5],fs,3),HbT_low,3);
    Behavior.R.gfp_HD_low_HbT_low = f_corr(f_bpf(gfp_HD,[0, 0.5],fs,3),HbT_low,3);
    Behavior.R.gfp_low_HbT_low = f_corr(f_bpf(gfp,[0, 0.5],fs,3),HbT_low,3);
    Behavior.signals = [squeeze(mean(rfp_HD.*comb_mask,[1, 2],'omitnan')),squeeze(mean(gfp_HD.*comb_mask,[1, 2],'omitnan')),squeeze(mean(HbT_low.*brain_mask,[1, 2],'omitnan')),Pupil,Whisking,Accelerometer];
    Behavior.R.signals = corrcoef(Behavior.signals);
    
    % NE IRF
    if string(mouseInfo.GRAB) == "GRAB_NE"
        [Behavior.NE_IRF.perf,Behavior.NE_IRF.IRF] = f_directDeco(f_bpf(gfp_HD,[0, 0.5],fs,3),rfp_HD,[-5, 10],fs,comb_mask,4);
    end
    
    clear gfp;
    
    %% Hb and HbO modeling
    
    if string(mouseInfo.GRAB) == "GRAB_NE"
        fprintf('\n\tHb modeling Analysis...');
        Hb_model = struct;
        
        ds = 32; % downsampling factor for estimation of initial timing parameters of IRF
        numThreads = 4; % number of cores to use in optimization algorithm
        win = [-5 10]; % IRF kernel range (s)
        
        Hb_model.irf_win = win;
        
        globalNE = mean(gfp_HD./std(gfp_HD,0,3).*brain_mask.*vessel_mask,[1,2],'omitnan');

        % HbO
        HbO_low = f_bpf(HbO,[0, 0.5],fs,3);
        [Hb_model.HbO.LR.perf,Hb_model.HbO.LR.params] = f_LR_varWeights(HbO_low,rfp_HD,globalNE.*ones(size(rfp_HD,[1,2])),win,fs,brain_mask,ds,numThreads);
        
        [Hb_model.HbO.IRFx2.perf,Hb_model.HbO.IRFx2.IRF,Hb_model.HbO.IRFx2.params] = f_2alphaDeconvolve(HbO_low,rfp_HD,globalNE.*ones(size(rfp_HD,[1,2])),win,fs,brain_mask,ds,numThreads);
        Hb_model.HbO.IRFx2.params.A = Hb_model.HbO.IRFx2.params.A*sum(Hb_model.HbO.IRFx2.IRF(:,1));
        Hb_model.HbO.IRFx2.params.B = Hb_model.HbO.IRFx2.params.B*sum(Hb_model.HbO.IRFx2.IRF(:,2));
        
        % Hb
        Hb_low = f_bpf(Hb,[0, 0.5],fs,3);
        [Hb_model.Hb.LR.perf,Hb_model.Hb.LR.params] = f_LR_varWeights(Hb_low,rfp_HD,globalNE.*ones(size(rfp_HD,[1,2])),win,fs,brain_mask,ds,numThreads);
        
        [Hb_model.Hb.IRFx2.perf,Hb_model.Hb.IRFx2.IRF,Hb_model.Hb.IRFx2.params] = f_2alphaDeconvolve(Hb_low,rfp_HD,globalNE.*ones(size(rfp_HD,[1,2])),win,fs,brain_mask,ds,numThreads);
        Hb_model.Hb.IRFx2.params.A = Hb_model.Hb.IRFx2.params.A*sum(Hb_model.Hb.IRFx2.IRF(:,1));
        Hb_model.Hb.IRFx2.params.B = Hb_model.Hb.IRFx2.params.B*sum(Hb_model.Hb.IRFx2.IRF(:,2));
        clear HbO_low Hb_low
    end
    
    %% unfiltered
    
    if string(mouseInfo.GRAB) == "GRAB_NE"
        fprintf('\n\tUnfiltered Analysis...');
        unfiltered = struct;
        
        ds = 32; % downsampling factor for estimation of initial timing parameters of IRF
        numThreads = 4; % number of cores to use in optimization algorithm
        win = [-5 10]; % IRF kernel range (s)
        
        unfiltered.irf_win = win;
        
        globalNE = mean(gfp_HD./std(gfp_HD,0,3).*brain_mask.*vessel_mask,[1,2],'omitnan');

        % estimate linear regression model
        [unfiltered.LR.perf,unfiltered.LR.params] = f_LR_varWeights(HbT,rfp_HD,globalNE.*ones(size(rfp_HD,[1,2])),win,fs,brain_mask,ds,numThreads);
        
        % estimate double IRF model - varying weights
        [unfiltered.IRFx2.perf,unfiltered.IRFx2.IRF,unfiltered.IRFx2.params] = f_2alphaDeconvolve(HbT,rfp_HD,globalNE.*ones(size(rfp_HD,[1,2])),win,fs,brain_mask,ds,numThreads);
        unfiltered.IRFx2.params.A = unfiltered.IRFx2.params.A*sum(unfiltered.IRFx2.IRF(:,1));
        unfiltered.IRFx2.params.B = unfiltered.IRFx2.params.B*sum(unfiltered.IRFx2.IRF(:,2));
    end
    
    %% shuffled NE
    
    if string(mouseInfo.GRAB) == "GRAB_NE"
        fprintf('\n\tShuffled Analysis...');
        shuffled = struct;
        
        ds = 32; % downsampling factor for estimation of initial timing parameters of IRF
        numThreads = 4; % number of cores to use in optimization algorithm
        win = [-5 10]; % IRF kernel range (s)
        
        shuffled.irf_win = win;
        
        shift = round(size(rfp_HD,3)/4);
        
        globalNE = mean(gfp_HD./std(gfp_HD,0,3).*brain_mask.*vessel_mask,[1,2],'omitnan');

        for i = 1:3
            % estimate linear regression model
            [shuffled.LR(i).perf,shuffled.LR(i).params] = f_LR_varWeights(HbT_low,rfp_HD,circshift(globalNE.*ones(size(rfp_HD,[1,2])),i*shift,3),win,fs,brain_mask,ds,numThreads);
            
            % estimate double IRF model - varying weights
            [shuffled.IRFx2(i).perf,shuffled.IRFx2(i).IRF,shuffled.IRFx2(i).params] = f_2alphaDeconvolve(HbT_low,rfp_HD,circshift(globalNE.*ones(size(rfp_HD,[1,2])),i*shift,3),win,fs,brain_mask,ds,numThreads);
            shuffled.IRFx2(i).params.A = shuffled.IRFx2(i).params.A*sum(shuffled.IRFx2(i).IRF(:,1));
            shuffled.IRFx2(i).params.B = shuffled.IRFx2(i).params.B*sum(shuffled.IRFx2(i).IRF(:,2));
        end
    
        % summarize shuffling
        shuffled.summary.LR.perf = mean(cat(3,shuffled.LR(:).perf),3);
        shuffled.summary.IRFx2.perf = mean(cat(3,shuffled.IRFx2(:).perf),3);
        shuffled.summary.IRFx2.IRF = mean(cat(3,shuffled.IRFx2(:).IRF),3);
    
        tmpParams = [shuffled.LR(:).params];
        shuffled.summary.LR.params.tA = mean([tmpParams(:).tA]);
        shuffled.summary.LR.params.tB = mean([tmpParams(:).tB]);
        shuffled.summary.LR.params.A = mean(cat(3,tmpParams(:).A),3);
        shuffled.summary.LR.params.B = mean(cat(3,tmpParams(:).B),3);
    
        tmpParams = [shuffled.IRFx2(:).params];
        shuffled.summary.IRFx2.params.A = mean(cat(3,tmpParams(:).A),3);
        shuffled.summary.IRFx2.params.B = mean(cat(3,tmpParams(:).B),3);
    end
    
    %% NE regression
    if string(mouseInfo.GRAB) == "GRAB_NE"
        fprintf('\n\tNE regression Analysis...');
        NE_reg = struct;
        
        globalNE = mean(gfp_HD./std(gfp_HD,0,3).*brain_mask.*vessel_mask,[1,2],'omitnan');

        HbT_reg = HbT_low./std(HbT_low,0,3);
        HbT_reg = HbT_reg-Fig2.LR.params.B.*globalNE.*ones(size(rfp_HD,[1,2]));
        NE_reg.HbT_reg = squeeze(mean(HbT_reg.*permute(allen_masks,[1 2 4 3]),[1, 2],'omitnan'));
        NE_reg.HbT_reg_global_LR = squeeze(globalNE) \ Fig3.HbT;
        NE_reg.HbT_reg_global = Fig3.HbT - NE_reg.HbT_reg_global_LR.*squeeze(globalNE);
        
        % connectivity
        win = Fig3.win;
        NE_reg.win = win;
    
        NE_reg.FC.HbT_reg = f_funConGram(NE_reg.HbT_reg,win*fs);
        NE_reg.lowNE_FC.HbT_reg = mean(NE_reg.FC.HbT_reg(:,:,lowNE),3);
        NE_reg.highNE_FC.HbT_reg = mean(NE_reg.FC.HbT_reg(:,:,highNE),3);

        NE_reg.FC.HbT_reg_global = f_funConGram(NE_reg.HbT_reg_global,win*fs);
        NE_reg.lowNE_FC.HbT_reg_global = mean(NE_reg.FC.HbT_reg_global(:,:,lowNE),3);
        NE_reg.highNE_FC.HbT_reg_global = mean(NE_reg.FC.HbT_reg_global(:,:,highNE),3);
        
        NE_reg.FC_r.Ca_HbT = f_corr(Fig3.FC.Ca,Fig3.FC.HbT,3);
        NE_reg.FC_r.Ca_HbT_reg = f_corr(Fig3.FC.Ca,NE_reg.FC.HbT_reg,3);
        NE_reg.FC_r.Ca_HbT_reg_global = f_corr(Fig3.FC.Ca,NE_reg.FC.HbT_reg_global,3);
    
        clear HbT_reg
    end
    
    %% frequency-dependent connectivity
    
    if string(mouseInfo.GRAB) == "GRAB_NE"
        fprintf('\n\tFrequency-dependent connectivity Analysis...');
        labels = {'low','medium','high'};
        fr_range = [0, 0.1;
                    0.1, 0.5;
                    0.5, 5];
        
        FC_fr = struct;
        
        for i = 1:3
            tmpFilt = f_bpf(rfp_HD,fr_range(i,:),fs,3);
            FC_fr.(labels{i}).Ca = squeeze(mean(tmpFilt./std(tmpFilt,0,3).*vessel_mask.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
            FC_fr.(labels{i}).FC.Ca = f_funConGram(FC_fr.(labels{i}).Ca,win*fs);
        
            FC_fr.(labels{i}).lowNE_FC.Ca = mean(FC_fr.(labels{i}).FC.Ca(:,:,lowNE),3);
            FC_fr.(labels{i}).highNE_FC.Ca = mean(FC_fr.(labels{i}).FC.Ca(:,:,highNE),3);
        
            if i ~= 3
                tmpFilt = f_bpf(HbT,fr_range(i,:),fs,3);
                FC_fr.(labels{i}).HbT = squeeze(mean(tmpFilt./std(tmpFilt,0,3).*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
                
                FC_fr.(labels{i}).FC.HbT = f_funConGram(FC_fr.(labels{i}).HbT,win*fs);
        
                FC_fr.(labels{i}).lowNE_FC.HbT = mean(FC_fr.(labels{i}).FC.HbT(:,:,lowNE),3);
                FC_fr.(labels{i}).highNE_FC.HbT = mean(FC_fr.(labels{i}).FC.HbT(:,:,highNE),3);
            end
        end
    end
    
    %% spectra
    
    if string(mouseInfo.GRAB) == "GRAB_NE"
        fprintf('\n\tSpectral Analysis...');
        spectra = struct;
        
        tmp = rfp_HD./std(rfp_HD,0,3);
        spectra.Ca_norm = squeeze(mean(tmp.*vessel_mask.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
        tmp = HbT./std(HbT,0,3);
        spectra.HbT_norm = squeeze(mean(tmp.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
        
        spectra.Ca = squeeze(mean(rfp_HD.*vessel_mask.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
        spectra.HbT = squeeze(mean(HbT.*permute(allen_masks,[1 2 4 3]),[1 2],'omitnan'));
        spectra.NE = squeeze(mean(gfp_HD.*vessel_mask.*brain_mask,[1, 2],'omitnan'));
        
        [spectra.SPG.Ca,spectra.SPG.t,spectra.SPG.f] = f_morlet(spectra.Ca_norm,fs,[0.01 fs/2],300);
        spectra.SPG.HbT = f_morlet(spectra.HbT_norm,fs,[0.01 fs/2],300);
    
        NE_bin = prctile(spectra.NE,[30 70]);
        low_idx = spectra.NE < NE_bin(1);
        high_idx = spectra.NE > NE_bin(2);
    
        spectra.low_NE.Ca = squeeze(mean(spectra.SPG.Ca(low_idx,:,:),1));
        spectra.low_NE.HbT = squeeze(mean(spectra.SPG.HbT(low_idx,:,:),1));
    
        spectra.high_NE.Ca = squeeze(mean(spectra.SPG.Ca(high_idx,:,:),1));
        spectra.high_NE.HbT = squeeze(mean(spectra.SPG.HbT(high_idx,:,:),1));
    
        [~,fIdx] = min(abs(spectra.SPG.f'-[0.1, 0.5]));
        
        spectra.low_NE.f_low.Ca = squeeze(mean(spectra.low_NE.Ca(1:fIdx(1),:)));
        spectra.low_NE.f_low.HbT = squeeze(mean(spectra.low_NE.HbT(1:fIdx(1),:)));
        spectra.low_NE.f_med.Ca = squeeze(mean(spectra.low_NE.Ca(fIdx(1):fIdx(2),:)));
        spectra.low_NE.f_med.HbT = squeeze(mean(spectra.low_NE.HbT(fIdx(1):fIdx(2),:)));
        spectra.low_NE.f_high.Ca = squeeze(mean(spectra.low_NE.Ca(fIdx(2):end,:)));
        
        spectra.high_NE.f_low.Ca = squeeze(mean(spectra.high_NE.Ca(1:fIdx(1),:)));
        spectra.high_NE.f_low.HbT = squeeze(mean(spectra.high_NE.HbT(1:fIdx(1),:)));
        spectra.high_NE.f_med.Ca = squeeze(mean(spectra.high_NE.Ca(fIdx(1):fIdx(2),:)));
        spectra.high_NE.f_med.HbT = squeeze(mean(spectra.high_NE.HbT(fIdx(1):fIdx(2),:)));
        spectra.high_NE.f_high.Ca = squeeze(mean(spectra.high_NE.Ca(fIdx(2):end,:)));
    
        clear tmp;
    end
    
    %% plot main figures

    plotFig1(Fig1,files.save);
    if string(mouseInfo.GRAB) == "GRAB_NE"
        plotFig2(Fig2,Fig1,shuffled,files.save);
        plotFig3(Fig3,files.save);
    end
    
    %% organize main figures
    
    if saveMetadata
        
        metadata = struct;
        metadata.Mouse = mouseInfo.ID;
        metadata.Date = sessionInfo.Date;
        metadata.Run = sessionInfo.Run;
        metadata.iRun = sessionInfo.interRun;
        metadata.settings.brain_mask = brain_mask;
        metadata.settings.vessel_mask = vessel_mask;
        metadata.settings.allen_masks = allen_masks;
        metadata.upd = '25-06-30';
        metadata.GRAB = mouseInfo.GRAB;
    
        savePath = strsplit(nwb_list(nwbI).Path,'/');
        savePath = char(savePath(end));
        savePath = savePath(1:end-4);
        savePath = [savePath '.mat'];
        savePath = fullfile(parentDir,'Analysis',savePath);
        
        if string(mouseInfo.GRAB) == "GRAB_NE"
            save(savePath,'metadata','Fig1','GRAB_FC','Fig2','Fig3','Behavior','Hb_model','unfiltered','shuffled','NE_reg','FC_fr','spectra','-v7.3');
        else
            save(savePath,'metadata','Fig1','GRAB_FC','Behavior','-v7.3');
        end
        fprintf('\n\tDone!\n');
    end
end
%% additional functions

function addPaths(mainDir)
    addDirs = genpath(mainDir);
    addDirs = strsplit(addDirs,pathsep);
    addDirs = addDirs(~cellfun(@(x) any(startsWith(strsplit(x,filesep), '.')), addDirs));
    addDirs = strjoin(addDirs,pathsep);
    addpath(addDirs);
end

function plotFig1(Fig1,savePath)
    
    parentDir = f_path();
    refBM = load(fullfile(parentDir,'Figures/plot_types/refAllen.mat'));
    refParcellation = refBM.refParcellation;
    refBM = refBM.refBM;

    plotBM = refBM;
    plotBM(:,1:300) = NaN;

    % B
    t = 0.1:0.1:size(Fig1.Ca_allen,1)/10;

    f = figure(Position=[100 100 1000 800]);
    tiledlayout(3,1);
    ax1 = nexttile;hold on;
    plot(t,Fig1.Ca_allen(:,5),color=c_Ca);
    plot(t,Fig1.GRAB_allen(:,5)-10,color=c_GRAB);
    plot(t,Fig1.HbT_allen(:,5)-20,color=c_HbT);
    legend('Ca^2^+',string(Fig1.GRAB_type),'HbT');
    ax1.XAxis.Visible = 'off';
    
    ax2 = nexttile;hold on;
    plot(t,Fig1.Ca_allen(:,3),color=c_Ca);
    plot(t,Fig1.GRAB_allen(:,3)-10,color=c_GRAB);
    plot(t,Fig1.HbT_allen(:,3)-20,color=c_HbT);
    ax2.XAxis.Visible = 'off';

    ax3 = nexttile;hold on;
    plot(t,Fig1.Pupil_full*0.4+0.6,color=c_pupil);
    plot(t,rescale(Fig1.Whisking,0.3,0.6),color=[0 0.7 0.7]);
    plot(t,rescale(Fig1.Accelerometer,0,0.3),color=[0 0 0]);
    legend('Pupil','Whisking','Movement');
    ax3.YAxis.Visible = 'off';
    xlabel('Time (s)');

    set(ax1,'FontSize',14);
    set(ax2,'FontSize',14);
    set(ax3,'FontSize',14);
    
    if savePath
        exportgraphics(f,fullfile(savePath,'Fig1B.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % C
    f = figure(Position=[100 100 900 400]);
    tiledlayout(1,2);
    nexttile;
    plot(0:0.1:10,Fig1.IRFx1_inv.IRF);
    box off;
    xlim([0 7]);
    xlabel('Time (s)');
    ylabel('a.u.');
    set(gca,'FontSize',14);

    nexttile;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig1.allen_masks,Fig1.IRFx1_inv.perf,0).*plotBM,cmp=cmpvir,bounds=[0 1],title='Global IRF',clabel='r');
    set(gca,'YDir','reverse');
    
    if savePath
        exportgraphics(f,fullfile(savePath,'Fig1C.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % D
    f = figure(Position=[100 100 900 400]);
    tiledlayout(1,2);
    nexttile;
    plot(0:0.1:10,Fig1.IRFx1_SSp.IRF);
    box off;
    xlim([0 7]);
    xlabel('Time (s)');
    ylabel('a.u.');
    set(gca,'FontSize',14);

    nexttile;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig1.allen_masks,Fig1.IRFx1_SSp.perf,0).*plotBM,cmp=cmpvir,bounds=[0 1],title='SSp IRF',clabel='r');
    set(gca,'YDir','reverse');
    
    if savePath
        exportgraphics(f,fullfile(savePath,'Fig1D.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % E
    f = figure(Position=[100 100 900 400]);
    tiledlayout(1,2);
    nexttile;hold on;
    plot(0:0.1:10,mean(Fig1.IRFx1_var.IRF_allen(:,4:5),2),color=c_Orange);
    plot(0:0.1:10,Fig1.IRFx1_var.IRF_allen(:,2),color=c_darkCyan);
    legend('SSp-tr/ll','MOs');
    box off;
    xlim([0 7]);
    xlabel('Time (s)');
    ylabel('a.u.');
    set(gca,'FontSize',14);

    nexttile;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig1.allen_masks,Fig1.IRFx1_var.perf,0).*plotBM,cmp=cmpvir,bounds=[0 1],title='Variant IRF',clabel='r');
    set(gca,'YDir','reverse');
    
    if savePath
        exportgraphics(f,fullfile(savePath,'Fig1E.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % F
    
    f = figure;
    bar(1:3,[mean(f_ImgReg_allen(refParcellation,Fig1.allen_masks,Fig1.IRFx1_inv.perf,0).*plotBM,'all','omitnan'),mean(f_ImgReg_allen(refParcellation,Fig1.allen_masks,Fig1.IRFx1_SSp.perf,0).*plotBM,'all','omitnan'),mean(f_ImgReg_allen(refParcellation,Fig1.allen_masks,Fig1.IRFx1_var.perf,0).*plotBM,'all','omitnan')],FaceColor=[1,1,1],EdgeColor=c_darkCyan,LineWidth=2);
    box off;
    ylabel('Model Performance (r)');
    set(gca,'XTickLabels',{'global','SSp','variant'},'FontSize',14);
    
    if savePath
        exportgraphics(f,fullfile(savePath,'Fig1F.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % G
    
    pS = struct;
    pS.xlim = [-4 4];

    lm = fitlm(mean(Fig1.GRAB,2),Fig1.IRFx1_SSp.perf_dt(:,2));
    lm = table2array(lm.Coefficients);

    f = figure;hold on;
    scatter(mean(Fig1.GRAB,2),Fig1.IRFx1_SSp.perf_dt(:,2),70,'filled',markerFaceColor=c_Orange);
    plot(pS.xlim,pS.xlim*lm(2,1)+lm(1,1),'-k','LineWidth',2);
    xlim(pS.xlim);
    ylim([-1 1]);
    xlabel('GRAB (\DeltaF/F)');
    ylabel('r_S_S_p_-_I_R_F(MOs)');
    set(gca,'FontSize',14);

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig1G.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % H
    
    f = figure;
    f_plotAllenMap(Fig1.SSp_perf_vs_GRAB,cmp=cmpbbr,cLabel='r',mask=plotBM,cRange=0.7*[-1,1]);
    
    if savePath
        exportgraphics(f,fullfile(savePath,'Fig1H.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % I - can't be calculated with only 1 subject/session

end

function plotFig2(Fig2,Fig1,shuffled,savePath)

    parentDir = f_path();
    refBM = load(fullfile(parentDir,'Figures/plot_types/refAllen.mat'));
    refParcellation = refBM.refParcellation;
    refBM = refBM.refBM;

    plotBM = refBM;
    plotBM(:,1:300) = NaN;

    % A
    
    t = 0.1:0.1:numel(Fig2.NE)/10;

    f = figure(Position=[100 100 1000 800]);
    tiledlayout(3,1);
    ax1 = nexttile;hold on;
    plot(t,Fig2.Ca_norm(:,12),color=c_Ca);
    plot(t,Fig2.NE-6,color=c_GRAB);
    plot(t,Fig2.HbT_norm(:,12)-12,color=[c_HbT 0.5]);
    ax1.XAxis.Visible = 'off';
    ylabel('Z-Score');

    ax2 = nexttile;hold on;
    plot(t,Fig2.LR_Ca(:,12),color=c_Ca);
    plot(t,Fig2.LR_NE(:,12)-6,color=c_GRAB);
    plot(t,Fig2.HbT_norm(:,12)-12,color=[c_HbT 0.5]);
    plot(t,Fig2.LR_Ca(:,12)+Fig2.LR_NE(:,12)-12,color=[0 0.8 0.8]);
    ax2.XAxis.Visible = 'off';
    ylabel('LR Model');

    ax3 = nexttile;hold on;
    plot(t,Fig2.IRFx2_Ca(:,12),color=c_Ca);
    plot(t,Fig2.IRFx2_NE(:,12)-6,color=c_GRAB);
    plot(t,Fig2.HbT_norm(:,12)-12,color=[c_HbT 0.5]);
    plot(t,Fig2.IRFx2_Ca(:,12)+Fig2.IRFx2_NE(:,12)-12,color=[0 0.8 0.8]);
    xlabel('Time (s)');
    ylabel('IRFx2 Model');

    set(ax1,'FontSize',14);
    set(ax2,'FontSize',14);
    set(ax3,'FontSize',14);
    
    if savePath
        exportgraphics(f,fullfile(savePath,'Fig2A.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % B
    f = figure(Position=[100 100 900 400]);
    tiledlayout(1,2);
    nexttile;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig2.LR.params.A,0).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='A');
    set(gca,'YDir','reverse');

    nexttile;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig2.LR.params.B,0).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='B');
    set(gca,'YDir','reverse');
    
    if savePath
        exportgraphics(f,fullfile(savePath,'Fig2B.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % C
    f = figure;hold on;
    bar(1,Fig2.LR.params.tA,FaceColor=[1,1,1],EdgeColor=c_Ca,LineWidth=2);
    bar(2,Fig2.LR.params.tB,FaceColor=[1,1,1],EdgeColor=c_GRAB,LineWidth=2);
    legend('t_A','t_B');
    ylabel('Lag (s)');
    set(gca,'FontSize',14);

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig2C.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % D

    f = figure;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig2.LR.perf,0).*plotBM,cmp=cmpvir,bounds=[0 1],title='LR Performance');
    set(gca,'YDir','reverse');

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig2D.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % E
    
    f = figure;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig2.LR.perf-Fig1.IRFx1_inv.perf,0).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR vs. Global IRF Performance');
    set(gca,'YDir','reverse');

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig2E.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % F
    f = figure(Position=[100 100 900 400]);
    tiledlayout(1,2);
    nexttile;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig2.IRFx2.params.A,0).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='A');
    set(gca,'YDir','reverse');

    nexttile;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig2.IRFx2.params.A,0).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='B');
    set(gca,'YDir','reverse');
    
    if savePath
        exportgraphics(f,fullfile(savePath,'Fig2F.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % G
    f = figure;hold on;
    plot(-5:0.1:10,Fig2.IRFx2.IRF(:,1),color=c_Ca);
    plot(-5:0.1:10,Fig2.IRFx2.IRF(:,2),color=c_GRAB);
    legend('IRF_C_a','IRF_N_E');
    xlabel('Time (s)');
    ylabel('a.u.');
    set(gca,'FontSize',14);

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig2G.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % H

    f = figure;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig2.IRFx2.perf,0).*plotBM,cmp=cmpvir,bounds=[0 1],title='LR Performance');
    set(gca,'YDir','reverse');

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig2H.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % I
    
    f = figure;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig2.IRFx2.perf-Fig1.IRFx1_inv.perf,0).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR vs. Global IRF Performance');
    set(gca,'YDir','reverse');

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig2I.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % J
    
    f = figure(Position=[100 100 900 400]);
    tiledlayout(1,2);
    nexttile;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig2.allen_masks,shuffled.summary.LR.perf,0).*plotBM,cmp=cmpvir,bounds=[0 1],title='Shuffled LR');
    set(gca,'YDir','reverse');

    nexttile;
    f_plotMap(f_ImgReg_allen(refParcellation,Fig2.allen_masks,shuffled.summary.IRFx2.perf,0).*plotBM,cmp=cmpvir,bounds=[0 1],title='Shuffled IRFx2');
    set(gca,'YDir','reverse');
    
    if savePath
        exportgraphics(f,fullfile(savePath,'Fig2J.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % K

    b = [];
    b(1) = mean(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig1.IRFx1_inv.perf,0).*plotBM,'all','omitnan');
    b(2) = mean(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig1.IRFx1_SSp.perf,0).*plotBM,'all','omitnan');
    b(3) = mean(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig1.IRFx1_var.perf,0).*plotBM,'all','omitnan');
    b(4) = mean(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig2.LR.perf,0).*plotBM,'all','omitnan');
    b(5) = mean(f_ImgReg_allen(refParcellation,Fig2.allen_masks,Fig2.IRFx2.perf,0).*plotBM,'all','omitnan');
    b(6) = mean(f_ImgReg_allen(refParcellation,Fig2.allen_masks,shuffled.summary.LR.perf,0).*plotBM,'all','omitnan');
    b(7) = mean(f_ImgReg_allen(refParcellation,Fig2.allen_masks,shuffled.summary.IRFx2.perf,0).*plotBM,'all','omitnan');

    f = figure;
    bar(1:7,b,FaceColor=[1,1,1],EdgeColor=c_darkCyan,LineWidth=2);
    box off;
    ylabel('Model Performance (r)');
    set(gca,'XTickLabels',{'global','SSp','variant','LR','IRFx2','shuffled LR','shuffled IRFx2'},'FontSize',14);

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig2K.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

end

function plotFig3(Fig3,savePath)

    parentDir = f_path();
    refBM = load(fullfile(parentDir,'Figures/plot_types/refAllen.mat'));
    refParcellation = refBM.refParcellation;
    refBM = refBM.refBM;

    plotBM = refBM;
    plotBM(:,1:300) = NaN;

    t = 0.1:0.1:size(Fig3.Ca,1)/10;

    % A
    f = figure(Position=[100 100 800 800]);
    tiledlayout(4,1);
    ax1 = nexttile;hold on;
    plot(t,Fig3.HbT(:,2),color=c_HbT*0.7);
    plot(t,Fig3.HbT(:,5)-10,color=c_HbT*1.2);
    ax1.XAxis.Visible = 'off';
    ylabel('HbT');
    
    ax2 = nexttile;hold on;
    plot(t,Fig3.Ca(:,2),color=c_Ca*0.7);
    plot(t,Fig3.Ca(:,5)-10,color=c_Ca*1.2);
    ax2.XAxis.Visible = 'off';
    ylabel('Ca^2^+');

    ax3 = nexttile;hold on;
    plot(t,Fig3.HbT(:,2),color=c_GRAB*0.7);
    plot(t,Fig3.HbT(:,5)-10,color=c_GRAB*1.2);
    ax3.XAxis.Visible = 'off';
    ylabel('NE');

    ax4 = nexttile;hold on;
    plot(t,Fig3.Pupil*0.4+0.6,color=c_pupil);
    plot(t,rescale(Fig3.Whisking,0.3,0.6),color=[0 0.7 0.7]);
    plot(t,rescale(Fig3.Accelerometer,0,0.3),color=[0 0 0]);
    legend('Pupil','Whisking','Movement');
    ax4.YAxis.Visible = 'off';
    xlabel('Time (s)');

    set(ax1,'FontSize',14);
    set(ax2,'FontSize',14);
    set(ax3,'FontSize',14);
    set(ax4,'FontSize',14);

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig3A.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % B
    
    f = figure(Position=[100 100 800 300]);hold on;
    plot(15:6:size(Fig3.Ca,1)/10-15,squeeze(Fig3.FC.Ca(2,5,:)),color=c_Ca);
    plot(15:6:size(Fig3.Ca,1)/10-15,squeeze(Fig3.FC.HbT(2,5,:)),color=c_HbT);
    ylim([0 1]);
    box off;
    xlabel('Time (s)');
    ylabel('r');
    legend('FC_C_a(MOs,SSp-ll)','FC_H_b_T(MOs,SSp-ll)')
    set(gca,'FontSize',14);
    
    if savePath
        exportgraphics(f,fullfile(savePath,'Fig3B.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % C

    f = figure(Position=[100 100 800 300]);
    plot(15:6:size(Fig3.Ca,1)/10-15,Fig3.FC.Ca_vs_HbT,color=c_darkCyan);
    ylim([0 1]);
    box off;
    xlabel('Time (s)');
    ylabel('r');
    set(gca,'FontSize',14);

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig3C.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % D
    
    pS = struct;
    pS.xlim = [-4 4];

    f = figure(Position=[100 100 900 400]);
    tiledlayout(1,2);
    nexttile;hold on;
    scatter(Fig3.GRAB,squeeze(Fig3.FC.Ca(2,5,:)),70,'filled',markerFaceColor=c_Orange);
    lm = fitlm(Fig3.GRAB,squeeze(Fig3.FC.Ca(2,5,:)));
    lm = table2array(lm.Coefficients);
    plot(pS.xlim,pS.xlim*lm(2,1)+lm(1,1),'-k','LineWidth',2);
    xlabel('Norepinephrine (\DeltaF/F)');
    ylabel('FC_C_a(MOs,SSp-ll)');
    set(gca,'FontSize',14);

    nexttile;hold on;
    scatter(Fig3.GRAB,squeeze(Fig3.FC.HbT(2,5,:)),70,'filled',markerFaceColor=c_Orange);
    lm = fitlm(Fig3.GRAB,squeeze(Fig3.FC.HbT(2,5,:)));
    lm = table2array(lm.Coefficients);
    plot(pS.xlim,pS.xlim*lm(2,1)+lm(1,1),'-k','LineWidth',2);
    xlabel('Norepinephrine (\DeltaF/F)');
    ylabel('FC_H_b_T(MOs,SSp-ll)');
    set(gca,'FontSize',14);
    
    if savePath
        exportgraphics(f,fullfile(savePath,'Fig3D.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % E
    
    pS = struct;
    pS.xlim = [-4 4];

    f = figure;hold on;
    scatter(Fig3.GRAB,Fig3.FC.Ca_vs_HbT,70,'filled',markerFaceColor=c_Orange);
    lm = fitlm(Fig3.GRAB,Fig3.FC.Ca_vs_HbT);
    lm = table2array(lm.Coefficients);
    plot(pS.xlim,pS.xlim*lm(2,1)+lm(1,1),'-k','LineWidth',2);
    xlabel('Norepinephrine (\DeltaF/F)');
    ylabel('r(FC_C_a,FC_H_b_T)');
    set(gca,'FontSize',14);

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig3E.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % F
    
    f = figure;hold on;
    bar(1,Fig3.r.FC_Ca_vs_GRAB(2,5),FaceColor=[1,1,1],EdgeColor=c_Ca,LineWidth=2);
    bar(2,Fig3.r.FC_HbT_vs_GRAB(2,5),FaceColor=[1,1,1],EdgeColor=c_HbT,LineWidth=2);
    bar(3,Fig3.r.FC_Ca_HbT_vs_GRAB,FaceColor=[1,1,1],EdgeColor=c_darkCyan,LineWidth=2);
    box off;
    ylabel('NE vs. FC (r)');
    legend('FC_C_a(MOs,SSp-ll)','FC_H_b_T(MOs,SSp-ll)','r(FC_C_a,FC_H_b_T)');
    set(gca,'FontSize',14);
    ax = gca;
    ax.XAxis.Visible = 'off';

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig3F.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % G
    f = figure(Position=[100 100 1300 400]);
    tiledlayout(1,3);
    nexttile;
    f_plotFC(Fig3.lowNE_FC.Ca,1,cmp=cmpvir,bounds=[0 1],title='Low NE Ca++ Connectivity',clabel='r');
    
    nexttile;
    f_plotFC(Fig3.highNE_FC.Ca,1,cmp=cmpvir,bounds=[0 1],title='High NE Ca++ Connectivity',clabel='r');
    
    ax3 = nexttile;
    f_plotFC(Fig3.highNE_FC.Ca-Fig3.lowNE_FC.Ca,1,cmp=cmpvir,bounds=0.25*[-1 1],title='High - Low NE Ca++ Connectivity',clabel='\Deltar');
    colormap(ax3,cmpbbr);

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig3G.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    % H
    f = figure(Position=[100 100 1300 400]);
    tiledlayout(1,3);
    nexttile;
    f_plotFC(Fig3.lowNE_FC.HbT,1,cmp=cmpvir,bounds=[0 1],title='Low NE HbT Connectivity',clabel='r');
    
    nexttile;
    f_plotFC(Fig3.highNE_FC.HbT,1,cmp=cmpvir,bounds=[0 1],title='High NE HbT Connectivity',clabel='r');
    
    ax3 = nexttile;
    f_plotFC(Fig3.highNE_FC.HbT-Fig3.lowNE_FC.HbT,1,cmp=cmpvir,bounds=0.25*[-1 1],title='High - Low NE HbT Connectivity',clabel='\Deltar');
    colormap(ax3,cmpbbr);

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig3H.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    %% I
    low = Fig3.lowNE_FC.Ca;
    high = Fig3.highNE_FC.Ca;
    
    combMap = [];
    tmpMask = NaN(size(refBM));
    
    combDiffMap = [];

    for i = [2,5,12]
        colMap = [];
        
        for rIdx = 1:12
            tmpMask(refParcellation.Masks(:,:,rIdx,2)) = low(i,rIdx);
        end
        colMap = tmpMask.*refBM;

        for rIdx = 1:12
            tmpMask(refParcellation.Masks(:,:,rIdx,2)) = high(i,rIdx);
        end
        colMap = [colMap;tmpMask.*refBM];
        colMap = colMap(:,301:end);

        combMap = [combMap colMap];

        for rIdx = 1:12
            tmpMask(refParcellation.Masks(:,:,rIdx,2)) = high(i,rIdx)-low(i,rIdx);
        end
        cropImg = tmpMask.*refBM;
        cropImg = cropImg(:,301:end);

        combDiffMap = [combDiffMap cropImg];

    end

    f = figure(Position=[100 100 600 800]);
    tiledlayout(3,1);
    ax1 = nexttile([2,1]);
    imagesc(combMap,AlphaData=~isnan(combMap));
    axis image off;
    c = colorbar;
    colormap cmpvir;
    c.Label.String = 'r';
    title('MOs        SSp-ll        VISp')
    clim([0 1]);

    ax2 = nexttile;
    imagesc(combDiffMap,AlphaData=~isnan(combDiffMap));
    axis image off;
    c = colorbar;
    clim(0.25*[-1 1]);
    c.Label.String = '\Deltar';
    colormap(ax2,cmpbbr);

    set(ax1,'FontSize',14);
    set(ax2,'FontSize',14);

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig3I.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

    %% J
    low = Fig3.lowNE_FC.HbT;
    high = Fig3.highNE_FC.HbT;
    
    combMap = [];
    tmpMask = NaN(size(refBM));
    
    combDiffMap = [];

    for i = [2,5,12]
        colMap = [];
        
        for rIdx = 1:12
            tmpMask(refParcellation.Masks(:,:,rIdx,2)) = low(i,rIdx);
        end
        colMap = tmpMask.*refBM;

        for rIdx = 1:12
            tmpMask(refParcellation.Masks(:,:,rIdx,2)) = high(i,rIdx);
        end
        colMap = [colMap;tmpMask.*refBM];
        colMap = colMap(:,301:end);

        combMap = [combMap colMap];

        for rIdx = 1:12
            tmpMask(refParcellation.Masks(:,:,rIdx,2)) = high(i,rIdx)-low(i,rIdx);
        end
        cropImg = tmpMask.*refBM;
        cropImg = cropImg(:,301:end);

        combDiffMap = [combDiffMap cropImg];

    end

    f = figure(Position=[100 100 600 800]);
    tiledlayout(3,1);
    ax1 = nexttile([2,1]);
    imagesc(combMap,AlphaData=~isnan(combMap));
    axis image off;
    c = colorbar;
    colormap cmpvir;
    c.Label.String = 'r';
    title('MOs        SSp-ll        VISp')
    clim([0 1]);

    ax2 = nexttile;
    imagesc(combDiffMap,AlphaData=~isnan(combDiffMap));
    axis image off;
    c = colorbar;
    clim(0.1*[-1 1]);
    c.Label.String = '\Deltar';
    colormap(ax2,cmpbbr);

    set(ax1,'FontSize',14);
    set(ax2,'FontSize',14);

    if savePath
        exportgraphics(f,fullfile(savePath,'Fig3J.png'),'Resolution',300,'BackgroundColor',[1 1 1]);
        close(f);
    end

end