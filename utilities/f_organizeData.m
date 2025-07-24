function [log, order, settings, Fig1, Fig2, Fig3, Behavior, Hb_model, unfiltered, shuffled, NE_reg, spectra, FC_fr, GRAB_FC] = f_organizeData(path)

files = dir(path);
files(1:2) = [];

log = struct('Mouse',[],'Date',[],'GRAB',[],'Run',[],'iRun',[]);

settings = struct;
Fig1 = struct;
Fig2 = struct;
Fig3 = struct;
Behavior = struct;
Hb_model = struct;
unfiltered = struct;
shuffled = struct;
NE_reg = struct;
spectra = struct;
FC_fr = struct;
GRAB_FC = struct;

N = numel(files);

for i = 1:N
    name = files(i).name;
    name = strsplit(name,'_');

    log(i).Mouse = strrep(name{1}(5:end),'-','_');
    log(i).Date = name{2}(5:end);
    log(i).Run = str2double(name{3}(5:6));
    log(i).iRun = str2double(name{4}(6:7));
    
    metadata = load(fullfile(files(i).folder,files(i).name));
    log(i).GRAB = metadata.metadata.GRAB;

    settings.brain_mask{i} = metadata.metadata.settings.brain_mask;
    settings.vessel_mask{i} = metadata.metadata.settings.vessel_mask;
    settings.allen_masks{i} = metadata.metadata.settings.allen_masks;
    settings.parcellation{i} = metadata.metadata.settings.parcellation.parcellation;

    Fig1.inv_perf{i} = metadata.Fig1.IRFx1_inv.perf;
    Fig1.inv_IRF{i} = metadata.Fig1.IRFx1_inv.IRF;
    Fig1.SSp_perf{i} = metadata.Fig1.IRFx1_SSp.perf;
    Fig1.SSp_IRF{i} = metadata.Fig1.IRFx1_SSp.IRF;
    Fig1.var_perf{i} = metadata.Fig1.IRFx1_var.perf;
    Fig1.var_IRF{i} = metadata.Fig1.IRFx1_var.IRF_allen;
    Fig1.GRAB{i} = metadata.Fig1.GRAB(:,2);
    Fig1.SSp_perf_dt{i} = metadata.Fig1.IRFx1_SSp.perf_dt(:,2);
    Fig1.SSp_perf_vs_GRAB{i} = metadata.Fig1.SSp_perf_vs_GRAB;
    
    corrWin = metadata.Fig1.perf_dt_win;
    tmpGRAB = movmean(metadata.GRAB_FC.GRAB_global,corrWin(2)*10,1);
    tmpGRAB = tmpGRAB(corrWin(1)*10/2:corrWin(2)*10:numel(tmpGRAB)-corrWin(1)*10/2);
    
    Fig1.GRAB_global{i} = tmpGRAB;
    Fig1.SSp_perf_vs_GRAB_global{i} = f_corr(metadata.Fig1.IRFx1_SSp.perf_dt,tmpGRAB,1);
    
    GRAB_FC.GRAB{i} = metadata.GRAB_FC.GRAB;
    GRAB_FC.GRAB_global{i} = metadata.GRAB_FC.GRAB_global;
    GRAB_FC.GRAB_norm{i} = metadata.GRAB_FC.GRAB_norm;
    GRAB_FC.FC{i} = metadata.GRAB_FC.FC;
    GRAB_FC.FC_detrend{i} = metadata.GRAB_FC.FC_detrend;

    Behavior.SPG.fr = metadata.Behavior.SPC.fr;
    Behavior.SPG.rfp_HD{i} = metadata.Behavior.SPC.rfp_HD;
    Behavior.SPG.gfp_HD{i} = metadata.Behavior.SPC.gfp_HD;
    Behavior.SPG.gfp{i} = metadata.Behavior.SPC.gfp;
    Behavior.SPG.HbT{i} = metadata.Behavior.SPC.HbT;
    Behavior.COH.rfp_HD_gfp_HD{i} = metadata.Behavior.COH.rfp_HD_gfp_HD;
    Behavior.COH.rfp_HD_HbT{i} = metadata.Behavior.COH.rfp_HD_HbT;
    Behavior.COH.gfp_HbT{i} = metadata.Behavior.COH.gfp_HbT;
    Behavior.COH.gfp_HD_HbT{i} = metadata.Behavior.COH.gfp_HD_HbT;
    Behavior.PHI.rfp_HD_gfp_HD{i} = metadata.Behavior.PHI.rfp_HD_gfp_HD;
    Behavior.PHI.rfp_HD_HbT{i} = metadata.Behavior.PHI.rfp_HD_HbT;
    Behavior.PHI.gfp_HbT{i} = metadata.Behavior.PHI.gfp_HbT;
    Behavior.PHI.gfp_HD_HbT{i} = metadata.Behavior.PHI.gfp_HD_HbT;
    Behavior.XC.rfp_HD_gfp_HD{i} = metadata.Behavior.XC.rfp_HD_gfp_HD;
    Behavior.XC.rfp_HD_HbT{i} = metadata.Behavior.XC.rfp_HD_HbT;
    Behavior.XC.gfp_HbT{i} = metadata.Behavior.XC.gfp_HbT;
    Behavior.XC.gfp_HD_HbT{i} = metadata.Behavior.XC.gfp_HD_HbT;
    Behavior.XC.lag = metadata.Behavior.XC.lag;
    Behavior.R.rfp_HD_low_HbT_low{i} = metadata.Behavior.R.rfp_HD_low_HbT_low;
    Behavior.R.gfp_HD_low_HbT_low{i} = metadata.Behavior.R.gfp_HD_low_HbT_low;
    Behavior.R.gfp_low_HbT_low{i} = metadata.Behavior.R.gfp_low_HbT_low;
    Behavior.R.rfp_HD_low_gfp_HD_low{i} = metadata.Behavior.R.rfp_HD_low_gfp_HD_low;
    Behavior.R.signals{i} = metadata.Behavior.R.signals;

    if string(log(i).GRAB) == "GRAB_NE"
        Behavior.NE_IRF_perf{i} = metadata.Behavior.NE_IRF.perf;
        Behavior.NE_IRF_IRF{i} = metadata.Behavior.NE_IRF.IRF;
        
        GRAB_FC.gfp_HD_vs_HbT_low{i} = [];
        tmpHbT = f_bpf(metadata.Fig3.HbT,[0,0.5],10);
        tmpNE = GRAB_FC.GRAB{i};
        for r = 1:12
            GRAB_FC.gfp_HD_vs_HbT_low{i}(:,r) = xcorr(tmpNE(:,r),tmpHbT(:,r),100,'normalized');
        end

        Fig2.LR_perf{i} = metadata.Fig2.LR.perf;
        Fig2.LR_tA{i} = metadata.Fig2.LR.params.tA;
        Fig2.LR_tB{i} = metadata.Fig2.LR.params.tB;
        Fig2.LR_A{i} = metadata.Fig2.LR.params.A;
        Fig2.LR_B{i} = metadata.Fig2.LR.params.B;
        Fig2.IRFx2_perf{i} = metadata.Fig2.IRFx2.perf;
        Fig2.IRFx2_IRF{i} = metadata.Fig2.IRFx2.IRF;
        Fig2.IRFx2_A{i} = metadata.Fig2.IRFx2.params.A;
        Fig2.IRFx2_B{i} = metadata.Fig2.IRFx2.params.B;
        Fig2.global_LR_perf{i} = metadata.Fig2.global.LR.perf;
        Fig2.global_LR_tA{i} = metadata.Fig2.global.LR.params.tA;
        Fig2.global_LR_tB{i} = metadata.Fig2.global.LR.params.tB;
        Fig2.global_LR_A{i} = metadata.Fig2.global.LR.params.A;
        Fig2.global_LR_B{i} = metadata.Fig2.global.LR.params.B;
        Fig2.global_IRFx2_perf{i} = metadata.Fig2.global.IRFx2.perf;
        Fig2.global_IRFx2_IRF{i} = metadata.Fig2.global.IRFx2.IRF;
        Fig2.global_IRFx2_A{i} = metadata.Fig2.global.IRFx2.params.A;
        Fig2.global_IRFx2_B{i} = metadata.Fig2.global.IRFx2.params.B;
        
        Fig3.NE{i} = metadata.Fig3.GRAB;
        Fig3.FC_Ca_MOs_SSpll{i} = squeeze(metadata.Fig3.FC.Ca(2,5,:));
        Fig3.FC_HbT_MOs_SSpll{i} = squeeze(metadata.Fig3.FC.HbT(2,5,:));
        Fig3.FC_Ca_vs_HbT{i} = metadata.Fig3.FC.Ca_vs_HbT;
        Fig3.lowNE_Ca{i} = metadata.Fig3.lowNE_FC.Ca;
        Fig3.lowNE_HbT{i} = metadata.Fig3.lowNE_FC.HbT;
        Fig3.highNE_Ca{i} = metadata.Fig3.highNE_FC.Ca;
        Fig3.highNE_HbT{i} = metadata.Fig3.highNE_FC.HbT;
        Fig3.FC_Ca_HbT_vs_GRAB{i} = metadata.Fig3.r.FC_Ca_HbT_vs_GRAB;
        Fig3.FC_Ca_vs_GRAB{i} = metadata.Fig3.r.FC_Ca_vs_GRAB;
        Fig3.FC_HbT_vs_GRAB{i} = metadata.Fig3.r.FC_HbT_vs_GRAB;

        Hb_model.Hb_LR_perf{i} = metadata.Hb_model.Hb.LR.perf;
        Hb_model.Hb_LR_tA{i} = metadata.Hb_model.Hb.LR.params.tA;
        Hb_model.Hb_LR_tB{i} = metadata.Hb_model.Hb.LR.params.tB;
        Hb_model.Hb_LR_A{i} = metadata.Hb_model.Hb.LR.params.A;
        Hb_model.Hb_LR_B{i} = metadata.Hb_model.Hb.LR.params.B;
        Hb_model.Hb_IRFx2_perf{i} = metadata.Hb_model.Hb.IRFx2.perf;
        Hb_model.Hb_IRFx2_IRF{i} = metadata.Hb_model.Hb.IRFx2.IRF;
        Hb_model.Hb_IRFx2_A{i} = metadata.Hb_model.Hb.IRFx2.params.A;
        Hb_model.Hb_IRFx2_B{i} = metadata.Hb_model.Hb.IRFx2.params.B;
        Hb_model.HbO_LR_perf{i} = metadata.Hb_model.HbO.LR.perf;
        Hb_model.HbO_LR_tA{i} = metadata.Hb_model.HbO.LR.params.tA;
        Hb_model.HbO_LR_tB{i} = metadata.Hb_model.HbO.LR.params.tB;
        Hb_model.HbO_LR_A{i} = metadata.Hb_model.HbO.LR.params.A;
        Hb_model.HbO_LR_B{i} = metadata.Hb_model.HbO.LR.params.B;
        Hb_model.HbO_IRFx2_perf{i} = metadata.Hb_model.HbO.IRFx2.perf;
        Hb_model.HbO_IRFx2_IRF{i} = metadata.Hb_model.HbO.IRFx2.IRF;
        Hb_model.HbO_IRFx2_A{i} = metadata.Hb_model.HbO.IRFx2.params.A;
        Hb_model.HbO_IRFx2_B{i} = metadata.Hb_model.HbO.IRFx2.params.B;

        unfiltered.LR_perf{i} = metadata.unfiltered.LR.perf;
        unfiltered.LR_tA{i} = metadata.unfiltered.LR.params.tA;
        unfiltered.LR_tB{i} = metadata.unfiltered.LR.params.tB;
        unfiltered.LR_A{i} = metadata.unfiltered.LR.params.A;
        unfiltered.LR_B{i} = metadata.unfiltered.LR.params.B;
        unfiltered.IRFx2_perf{i} = metadata.unfiltered.IRFx2.perf;
        unfiltered.IRFx2_IRF{i} = metadata.unfiltered.IRFx2.IRF;
        unfiltered.IRFx2_A{i} = metadata.unfiltered.IRFx2.params.A;
        unfiltered.IRFx2_B{i} = metadata.unfiltered.IRFx2.params.B;

        shuffled.LR_perf{i} = metadata.shuffled.summary.LR.perf;
        shuffled.LR_tA{i} = metadata.shuffled.summary.LR.params.tA;
        shuffled.LR_tB{i} = metadata.shuffled.summary.LR.params.tB;
        shuffled.LR_A{i} = metadata.shuffled.summary.LR.params.A;
        shuffled.LR_B{i} = metadata.shuffled.summary.LR.params.B;
        shuffled.IRFx2_perf{i} = metadata.shuffled.summary.IRFx2.perf;
        shuffled.IRFx2_IRF{i} = metadata.shuffled.summary.IRFx2.IRF;
        shuffled.IRFx2_A{i} = metadata.shuffled.summary.IRFx2.params.A;
        shuffled.IRFx2_B{i} = metadata.shuffled.summary.IRFx2.params.B;

        NE_reg.lowNE_FC{i} = metadata.NE_reg.lowNE_FC.HbT_reg;
        NE_reg.highNE_FC{i} = metadata.NE_reg.highNE_FC.HbT_reg;
        NE_reg.FC_R_HbT{i} = metadata.NE_reg.FC_r.Ca_HbT;
        NE_reg.FC_R_HbT_reg{i} = metadata.NE_reg.FC_r.Ca_HbT_reg;
        NE_reg.IRFx1_inv_perf{i} = metadata.NE_reg.IRFx1_inv.perf;
        NE_reg.IRFx1_inv_IRF{i} = metadata.NE_reg.IRFx1_inv.IRF;
        
        bounds = prctile(metadata.Fig3.GRAB,[30, 70]);
        lowNE = metadata.Fig3.GRAB < bounds(1);
        highNE = metadata.Fig3.GRAB > bounds(2);
        
        global_HbT_reg = metadata.GRAB_FC.GRAB_global \ metadata.Fig3.HbT;
        global_HbT_reg = metadata.Fig3.HbT-global_HbT_reg.*metadata.GRAB_FC.GRAB_global;

        global_FC_HbT_reg = f_funConGram(global_HbT_reg,metadata.Fig3.win*10);
        NE_reg.global_lowNE_FC{i} = mean(global_FC_HbT_reg(:,:,lowNE),3);
        NE_reg.global_highNE_FC{i} = mean(global_FC_HbT_reg(:,:,highNE),3);
        
        NE_reg.global_FC_R_HbT_reg{i} = f_corr(metadata.Fig3.FC.Ca,global_FC_HbT_reg,3);
        
        FC_fr.lowF_lowNE_Ca{i} = metadata.FC_fr.low.lowNE_FC.Ca;
        FC_fr.lowF_highNE_Ca{i} = metadata.FC_fr.low.highNE_FC.Ca;
        FC_fr.medF_lowNE_Ca{i} = metadata.FC_fr.medium.lowNE_FC.Ca;
        FC_fr.medF_highNE_Ca{i} = metadata.FC_fr.medium.highNE_FC.Ca;
        FC_fr.highF_lowNE_Ca{i} = metadata.FC_fr.high.lowNE_FC.Ca;
        FC_fr.highF_highNE_Ca{i} = metadata.FC_fr.high.highNE_FC.Ca;
        FC_fr.lowF_lowNE_HbT{i} = metadata.FC_fr.low.lowNE_FC.HbT;
        FC_fr.lowF_highNE_HbT{i} = metadata.FC_fr.low.highNE_FC.HbT;
        FC_fr.medF_lowNE_HbT{i} = metadata.FC_fr.medium.lowNE_FC.HbT;
        FC_fr.medF_highNE_HbT{i} = metadata.FC_fr.medium.highNE_FC.HbT;
        
        NE_bin = prctile(metadata.spectra.NE,[30 70]);
        low_idx = metadata.spectra.NE < NE_bin(1);
        high_idx = metadata.spectra.NE > NE_bin(2);
        spectra.fr = metadata.spectra.SPG.f';
        SPG_Ca = metadata.spectra.SPG.Ca.*metadata.spectra.SPG.f;
        SPG_HbT = metadata.spectra.SPG.HbT.*metadata.spectra.SPG.f;
        
        [~,fIdx] = min(abs(metadata.spectra.SPG.f'-[0.1, 0.5]));

        spectra.low_NE_Ca{i} = mean(metadata.spectra.low_NE.Ca,2);
        spectra.high_NE_Ca{i} = mean(metadata.spectra.high_NE.Ca,2);
        spectra.low_NE_HbT{i} = mean(metadata.spectra.low_NE.HbT,2);
        spectra.high_NE_HbT{i} = mean(metadata.spectra.high_NE.HbT,2);
        spectra.NE{i} = metadata.spectra.NE;
        spectra.SPG_Ca{i} = SPG_Ca;
        spectra.SPG_HbT{i} = SPG_HbT;

        spectra.lowF_lowNE_Ca{i} = squeeze(mean(SPG_Ca(low_idx,1:fIdx(1),:),[1,2]));
        spectra.lowF_highNE_Ca{i} = squeeze(mean(SPG_Ca(high_idx,1:fIdx(1),:),[1,2]));
        spectra.medF_lowNE_Ca{i} = squeeze(mean(SPG_Ca(low_idx,fIdx(1):fIdx(2),:),[1,2]));
        spectra.medF_highNE_Ca{i} = squeeze(mean(SPG_Ca(high_idx,fIdx(1):fIdx(2),:),[1,2]));
        spectra.highF_lowNE_Ca{i} = squeeze(mean(SPG_Ca(low_idx,fIdx(2):end,:),[1,2]));
        spectra.highF_highNE_Ca{i} = squeeze(mean(SPG_Ca(high_idx,fIdx(2):end,:),[1,2]));
        spectra.lowF_lowNE_HbT{i} = squeeze(mean(SPG_HbT(low_idx,1:fIdx(1),:),[1,2]));
        spectra.lowF_highNE_HbT{i} = squeeze(mean(SPG_HbT(high_idx,1:fIdx(1),:),[1,2]));
        spectra.medF_lowNE_HbT{i} = squeeze(mean(SPG_HbT(low_idx,fIdx(1):fIdx(2),:),[1,2]));
        spectra.medF_highNE_HbT{i} = squeeze(mean(SPG_HbT(high_idx,fIdx(1):fIdx(2),:),[1,2]));
        
    end

end

mice = unique({log.Mouse})';
M = numel(mice);

order = struct('Mouse',[],'Runs',[],'GRAB',[],'N',[]);

for mIdx = 1:M
    order(mIdx).Mouse = mice{mIdx};
    order(mIdx).Runs = find(strcmp({log.Mouse},mice{mIdx}));
    order(mIdx).GRAB = log(order(mIdx).Runs(1)).GRAB;
    order(mIdx).N = numel(order(mIdx).Runs);
end

end