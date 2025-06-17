function [perf,IRF,outParams] = f_2alphaDeconvolve(HbT,Ca,NE,win,fs, ...
    brain_mask,ds,maxThreads,initialParam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs deconvolution of HbT by applying a single alpha function 
% convolution kernel to Ca and NE, individually. Allows for spatially
% variant weights (A, B) for both alpha functions and spatially invariant
% timing parameters (tA,tB,tauA,tauB).
% INPUTS:
%   HbT - HbT video
%   Ca - Ca video
%   NE - NE video
%   win - window to use for IRF kernel (s) [t1 t2]
%   fs - sampling frequency (Hz)
%   brain_mask - mask of brain exposure (2D NaN image)
%   ds - downsampling kernel (int, typically 32)
%   maxThreads - maximum number of cores to use (default 4)
%   initialParams - initial timing parameters to use (default defined 
%       below)
% OUTPUTS:
%   perf - correlation map showing correlation between input HbT and
%       predicted HbT
%   IRF - both estimated IRFs
%   outParams - output parameters including tA, tB, tauA, tauB and brain
%       maps of A and B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 7
    maxNumCompThreads(maxThreads);
else
    maxNumCompThreads(4);
end

%% downsample and reshape data
dsHbT = f_downsample(HbT,ds);
dsCa = f_downsample(Ca,ds);
dsNE = f_downsample(NE,ds);

ds_brain_mask = f_downsample(brain_mask,ds);
nanIdx = ~isnan(ds_brain_mask(:));
N = sum(nanIdx);

dim = size(dsHbT);
dsHbT = reshape(dsHbT,[],dim(3));
dsCa = reshape(dsCa,[],dim(3));
dsNE = reshape(dsNE,[],dim(3));

dsHbT = dsHbT(nanIdx,:)';
dsCa = dsCa(nanIdx,:)';
dsNE = dsNE(nanIdx,:)';

dsHbT = dsHbT./std(dsHbT,0,1);
dsCa = dsCa./std(dsCa,0,1);
dsNE = dsNE./std(dsNE,0,1);

%% create design matrices
T = size(dsCa,1);
l_irf = fs*range(win)+1;
idx_irf = win(1)*fs:win(2)*fs;

i1 = abs(min([idx_irf;zeros(1,numel(idx_irf))],[],1))+1;
i2 = T-max([idx_irf;zeros(1,numel(idx_irf))],[],1);
i3 = max([idx_irf;zeros(1,numel(idx_irf))],[],1)+1;
i4 = T-i1+1;

Ca_mat = zeros(T,N,l_irf);
NE_mat = zeros(T,N,l_irf);

for v = 1:l_irf
    Ca_mat(i3(v):i4(v),:,v) = dsCa(i1(v):i2(v),:);
end
for v = 1:l_irf
    NE_mat(i3(v):i4(v),:,v) = dsNE(i1(v):i2(v),:);
end

Ca_mat = permute(Ca_mat,[3,2,1]);
NE_mat = permute(NE_mat,[3,2,1]);

design_matrix = [Ca_mat;NE_mat];
design_HbT = permute(dsHbT,[3 2 1]);

%% optimize initial parameters

if nargin > 8
    t1 = initialParam(1);
    t2 = initialParam(2);
    tau1 = initialParam(3);
    tau2 = initialParam(4);
else
    t1 = 0.1;
    t2 = -1;
    tau1 = 0.5;
    tau2 = 0.53;
end
A = ones(1,N);
B = ones(1,N);

hrf1 = f_alpha_IRF(t1,tau1,0.5,A,0,fs,win);
hrf2 = f_alpha_IRF(t2,tau2,0.5,B,0,fs,win);

conv_Ca = sum(Ca_mat.*hrf1,1);
conv_NE = sum(NE_mat.*hrf2,1);

LR = f_hemRegress(design_HbT,cat(4,conv_Ca,conv_NE),ones(1,N));

A = LR(:,:,1);
B = LR(:,:,2);

%% run optimization using gradient descent

fun = @(params)f_hrf_cost_func(params(1),params(2),params(3),params(4),params(5:N+4),params(N+5:2*N+4),fs,win,design_matrix,design_HbT);

options = optimset('MaxFunEvals',25000,'MaxIter',500,'Display','iter','Algorithm','active-set');
params = fmincon(fun,[t1, t2, tau1, tau2, A, B],[],[],[],[],[0 win(1) 0.01 0.01 -Inf*ones(1,numel(A)*2)],[win(2) win(2) win(2) win(2) Inf*ones(1,numel(A)*2)],[],options);

%% pixel-wise LR

hrf1 = f_alpha_IRF(params(1),params(3),0.5,1,0,fs,win);
hrf2 = f_alpha_IRF(params(2),params(4),0.5,1,0,fs,win);

conv_Ca = f_3Dconvolve(Ca./std(Ca,0,3),hrf1,win*fs,brain_mask);
conv_NE = f_3Dconvolve(NE./std(NE,0,3),hrf2,win*fs,brain_mask);

LR = f_hemRegress(HbT./std(HbT,0,3),cat(4,conv_Ca,conv_NE),brain_mask);

outParams = struct;
outParams.tA = params(1);
outParams.tB = params(2);
outParams.tauA = params(3);
outParams.tauB = params(4);
outParams.A = LR(:,:,1);
outParams.B = LR(:,:,2);

IRF = [hrf1,-hrf2];

%% calculate performance
pred_HbT = conv_Ca.*LR(:,:,1)+conv_NE.*LR(:,:,2);

perf = f_corr(HbT./std(HbT,0,3),pred_HbT,3);

%%

function J = f_hrf_cost_func(t1,t2,tau1,tau2,A,B,sr,range,X_mat,y)
    HRF1 = f_alpha_IRF(t1,tau1,0.5,A,0,sr,range);
    HRF2 = f_alpha_IRF(t2,tau2,0.5,B,0,sr,range);

    hrf = [HRF1; HRF2];
    conv_result = sum(X_mat.*hrf,1);
    J = sqrt(mean((y - conv_result).^2,'all'));
end

function IRF = f_alpha_IRF(t0,tau1,tau2,A,B,sr,range)
    hrf_l = range*sr;
    tr = (((hrf_l(1):hrf_l(2))/sr)-t0)';
    D = ((tr)./tau1).^3.*exp(-(tr)./tau1);
    D(tr<0) = 0;
    C = ((tr)./tau2).^3.*exp(-(tr)./tau2);
    C(tr<0) = 0;
    
    IRF = A.*D + B.*C;
end

end