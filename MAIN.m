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

%% check for analyzed runs

N = numel(nwb_list);
processed = false(N,1);

for i = 1:N
    savePath = strsplit(nwb_list(i).Path,'/');
    savePath = char(savePath(end));
    savePath = savePath(1:end-4);
    savePath = [savePath '.mat'];
    savePath = fullfile('/projectnb/devorlab/bcraus/AnalysisCode/Neuromodulation/Analysis',savePath);

    if f_checkFile(savePath)
        processed(i) = true;
    end
end

indices = find(~processed)';

%%
% temporary index to analyze only a single run
for nwbI = indices

    fprintf(['Analyzing ', nwb_list(nwbI).Path, '...']);

    nwb = nwbRead(nwb_list(nwbI).Path);
    
    [~,gfp,rfp_HD,gfp_HD,Hb,HbO,HbT,Whisking,Pupil,Accelerometer,brain_mask,vessel_mask,allen_masks,fs,mouseInfo,sessionInfo] = f_extractNWB(nwb);
    
    GRAB = mouseInfo.GRAB;
    
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
    HbT_low = f_bpf(HbT,[0, 0.5],fs,3);
    
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
    tmp = detrend(GRAB_FC.GRAB);
    GRAB_FC.FC_detrend = corrcoef(tmp);

    %% perform analysis for Fig 2
    
    if string(GRAB) == "GRAB_NE"
        fprintf('\n\tFigure 2 Analysis...');
        Fig2 = struct;
        
        ds = 32; % downsampling factor for estimation of initial timing parameters of IRF
        % numThreads = 4; % number of cores to use in optimization algorithm
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

        % repeat for global NE
        
        globalNE = mean(gfp_HD./std(gfp_HD,0,3).*brain_mask.*vessel_mask,[1,2],'omitnan');

        Fig2.global.LR = struct;
        Fig2.global.NE = squeeze(globalNE);

        [Fig2.global.LR.perf,Fig2.global.LR.params] = f_LR_varWeights(HbT_low,rfp_HD,globalNE.*ones(size(rfp_HD,[1,2])),win,fs,brain_mask,ds,numThreads);
        
        % estimate double IRF model - varying weights
        
        Fig2.global.IRFx2 = struct;
        [Fig2.global.IRFx2.perf,Fig2.global.IRFx2.IRF,Fig2.global.IRFx2.params] = f_2alphaDeconvolve(HbT_low,rfp_HD,globalNE.*ones(size(rfp_HD,[1,2])),win,fs,brain_mask,ds,numThreads);
        Fig2.global.IRFx2.params.A = Fig2.global.IRFx2.params.A*sum(Fig2.global.IRFx2.IRF(:,1));
        Fig2.global.IRFx2.params.B = Fig2.global.IRFx2.params.B*sum(Fig2.global.IRFx2.IRF(:,2));
    
    end
    
    %% perform analysis for Fig 3
    
    if string(GRAB) == "GRAB_NE"
        fprintf('\n\tFigure 3 Analysis...');
        Fig3 = struct;
        
        win = [30 6]; % length of sliding connectivity window (s)
        Fig3.win = win;
        
        % extract Ca and HbT ROIs
        
        HbT_low = f_bpf(HbT,[0, 0.5],fs,3);
        
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
    if string(GRAB) == "GRAB_NE"
        [Behavior.NE_IRF.perf,Behavior.NE_IRF.IRF] = f_directDeco(f_bpf(gfp_HD,[0, 0.5],fs,3),rfp_HD,[-5, 10],fs,comb_mask,4);
    end
    
    clear gfp;
    
    %% Hb and HbO modeling
    
    if string(GRAB) == "GRAB_NE"
        fprintf('\n\tHb modeling Analysis...');
        Hb_model = struct;
        
        ds = 32; % downsampling factor for estimation of initial timing parameters of IRF
        numThreads = 4; % number of cores to use in optimization algorithm
        win = [-5 10]; % IRF kernel range (s)
        
        Hb_model.irf_win = win;
        
        % HbO
        HbO_low = f_bpf(HbO,[0, 0.5],fs,3);
        [Hb_model.HbO.LR.perf,Hb_model.HbO.LR.params] = f_LR_varWeights(HbO_low,rfp_HD,gfp_HD,win,fs,brain_mask,ds,numThreads);
        
        [Hb_model.HbO.IRFx2.perf,Hb_model.HbO.IRFx2.IRF,Hb_model.HbO.IRFx2.params] = f_2alphaDeconvolve(HbO_low,rfp_HD,gfp_HD,win,fs,brain_mask,ds,numThreads);
        Hb_model.HbO.IRFx2.params.A = Hb_model.HbO.IRFx2.params.A*sum(Hb_model.HbO.IRFx2.IRF(:,1));
        Hb_model.HbO.IRFx2.params.B = Hb_model.HbO.IRFx2.params.B*sum(Hb_model.HbO.IRFx2.IRF(:,2));
        
        % Hb
        Hb_low = f_bpf(Hb,[0, 0.5],fs,3);
        [Hb_model.Hb.LR.perf,Hb_model.Hb.LR.params] = f_LR_varWeights(Hb_low,rfp_HD,gfp_HD,win,fs,brain_mask,ds,numThreads);
        
        [Hb_model.Hb.IRFx2.perf,Hb_model.Hb.IRFx2.IRF,Hb_model.Hb.IRFx2.params] = f_2alphaDeconvolve(Hb_low,rfp_HD,gfp_HD,win,fs,brain_mask,ds,numThreads);
        Hb_model.Hb.IRFx2.params.A = Hb_model.Hb.IRFx2.params.A*sum(Hb_model.Hb.IRFx2.IRF(:,1));
        Hb_model.Hb.IRFx2.params.B = Hb_model.Hb.IRFx2.params.B*sum(Hb_model.Hb.IRFx2.IRF(:,2));
        clear HbO_low Hb_low
    end
    
    %% unfiltered
    
    if string(GRAB) == "GRAB_NE"
        fprintf('\n\tUnfiltered Analysis...');
        unfiltered = struct;
        
        ds = 32; % downsampling factor for estimation of initial timing parameters of IRF
        numThreads = 4; % number of cores to use in optimization algorithm
        win = [-5 10]; % IRF kernel range (s)
        
        unfiltered.irf_win = win;
        
        % estimate linear regression model
        [unfiltered.LR.perf,unfiltered.LR.params] = f_LR_varWeights(HbT,rfp_HD,gfp_HD,win,fs,brain_mask,ds,numThreads);
        
        % estimate double IRF model - varying weights
        [unfiltered.IRFx2.perf,unfiltered.IRFx2.IRF,unfiltered.IRFx2.params] = f_2alphaDeconvolve(HbT,rfp_HD,gfp_HD,win,fs,brain_mask,ds,numThreads);
        unfiltered.IRFx2.params.A = unfiltered.IRFx2.params.A*sum(unfiltered.IRFx2.IRF(:,1));
        unfiltered.IRFx2.params.B = unfiltered.IRFx2.params.B*sum(unfiltered.IRFx2.IRF(:,2));
    end
    
    %% shuffled NE
    
    if string(GRAB) == "GRAB_NE"
        fprintf('\n\tShuffled Analysis...');
        shuffled = struct;
        
        ds = 32; % downsampling factor for estimation of initial timing parameters of IRF
        numThreads = 4; % number of cores to use in optimization algorithm
        win = [-5 10]; % IRF kernel range (s)
        
        shuffled.irf_win = win;
            
        shift = round(size(rfp_HD,3)/4);
        
        for i = 1:3
            % estimate linear regression model
            [shuffled.LR(i).perf,shuffled.LR(i).params] = f_LR_varWeights(HbT_low,rfp_HD,circshift(gfp_HD,i*shift,3),win,fs,brain_mask,ds,numThreads);
            
            % estimate double IRF model - varying weights
            [shuffled.IRFx2(i).perf,shuffled.IRFx2(i).IRF,shuffled.IRFx2(i).params] = f_2alphaDeconvolve(HbT_low,rfp_HD,circshift(gfp_HD,i*shift,3),win,fs,brain_mask,ds,numThreads);
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
    if string(GRAB) == "GRAB_NE"
        fprintf('\n\tNE regression Analysis...');
        NE_reg = struct;
        
        HbT_reg = HbT_low./std(HbT_low,0,3);
        HbT_reg = HbT_reg-Fig2.LR.params.B.*gfp_HD./std(gfp_HD,0,3);
        NE_reg.HbT_reg = squeeze(mean(HbT_reg.*permute(allen_masks,[1 2 4 3]),[1, 2],'omitnan'));
        
        [NE_reg.IRFx1_inv.perf,NE_reg.IRFx1_inv.IRF,NE_reg.IRFx1_inv.params,tmpCorr] = f_1xIRF(HbT_reg,rfp_HD,win,fs,brain_mask,4,corrWin*fs,numThreads);
        
        % connectivity
        win = Fig3.win;
        NE_reg.win = win;
    
        NE_reg.FC.HbT_reg = f_funConGram(NE_reg.HbT_reg,win*fs);
        NE_reg.lowNE_FC.HbT_reg = mean(NE_reg.FC.HbT_reg(:,:,lowNE),3);
        NE_reg.highNE_FC.HbT_reg = mean(NE_reg.FC.HbT_reg(:,:,highNE),3);
        
        NE_reg.FC_r.Ca_HbT = f_corr(Fig3.FC.Ca,Fig3.FC.HbT,3);
        NE_reg.FC_r.Ca_HbT_reg = f_corr(Fig3.FC.Ca,NE_reg.FC.HbT_reg,3);
    
        clear HbT_reg
    end
    
    %% frequency-dependent connectivity
    
    if string(GRAB) == "GRAB_NE"
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
    
    if string(GRAB) == "GRAB_NE"
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
    
    %% organize main figures
    
    parcellation = load(fullfile('/projectnb/devorlab/bcraus/HRF/1P',sessionInfo.Date,mouseInfo.ID,'DataAnalysis','parcells_upd'));
    
    metadata = struct;
    metadata.Mouse = mouseInfo.ID;
    metadata.Date = sessionInfo.Date;
    metadata.Run = sessionInfo.Run;
    metadata.iRun = sessionInfo.interRun;
    metadata.settings.brain_mask = brain_mask;
    metadata.settings.vessel_mask = vessel_mask;
    metadata.settings.allen_masks = allen_masks;
    metadata.settings.parcellation = parcellation;
    metadata.upd = '25-06-30';
    metadata.GRAB = GRAB;

    savePath = strsplit(nwb_list(nwbI).Path,'/');
    savePath = char(savePath(end));
    savePath = savePath(1:end-4);
    savePath = [savePath '.mat'];
    savePath = fullfile('/projectnb/devorlab/bcraus/AnalysisCode/Neuromodulation/Analysis',savePath);
    
    if string(GRAB) == "GRAB_NE"
        save(savePath,'metadata','Fig1','Fig2','Fig3','Behavior','Hb_model','unfiltered','shuffled','NE_reg','FC_fr','spectra','-v7.3');
    else
        save(savePath,'metadata','Fig1','Behavior','-v7.3');
    end
    fprintf('\n\tDone!\n');

end
%% additional functions

function addPaths(mainDir)
    addDirs = genpath(mainDir);
    addDirs = strsplit(addDirs,pathsep);
    addDirs = addDirs(~cellfun(@(x) any(startsWith(strsplit(x,filesep), '.')), addDirs));
    addDirs = strjoin(addDirs,pathsep);
    addpath(addDirs);
end
