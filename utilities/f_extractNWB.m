function [rfp_HD,gfp_HD,Hb,HbO,HbT,Whisking,Pupil,Accelerometer,brain_mask,allen_masks,fs] = f_extractNWB(nwb)

%% extract optical channels

rotation = 270; % convert between python and MATLAB

rfp_HD = 100*imrotate(nwb.acquisition.get('rfp_HD').data.load(),rotation);
gfp_HD = 100*imrotate(nwb.acquisition.get('gfp_HD').data.load(),rotation);
Hb = 1e6*imrotate(nwb.acquisition.get('Hb').data.load(),rotation);
HbO = 1e6*imrotate(nwb.acquisition.get('HbO').data.load(),rotation);
HbT = Hb+HbO;

%% extract behavioral readouts

Whisking = nwb.processing.get('behavior').nwbdatainterface.get('Behavior').timeseries.get('Whisking').data.load();
Pupil = nwb.processing.get('behavior').nwbdatainterface.get('Behavior').timeseries.get('Pupil diameter').data.load();
Accelerometer = nwb.processing.get('behavior').nwbdatainterface.get('Behavior').timeseries.get('Accelerometer').data.load();

%% extract masks

brain_mask = imrotate(nwb.processing.get('ophys').nwbdatainterface.get('ImageSegmentation').planesegmentation.get('brain_mask').image_mask.data.load(),rotation);
allen_masks = double(imrotate(nwb.processing.get('ophys').nwbdatainterface.get('ImageSegmentation').planesegmentation.get('allen_masks').image_mask.data.load(),rotation));
allen_masks(allen_masks==0) = NaN;

%% extract framerate

fs = nwb.acquisition.get('Hb').starting_time_rate;

end