function [perf,IRF,outParams,corrGram,pred_HbT] = f_1xIRF_varWeights(HbT,Ca,win,fs, ...
    brain_mask,ds,corrWin,maxThreads,initialParam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs deconvolution of HbT by applying a single double alpha function 
% convolution kernel to Ca. Allows for spatially variant weights (A, B) for
% both alpha functions and spatially invariant timing parameters 
% (t0,tauA,tauB).
% INPUTS:
%   HbT - HbT video
%   Ca - Ca video
%   win - window to use for IRF kernel (s) [t1 t2]
%   fs - sampling frequency (Hz)
%   brain_mask - mask of brain exposure (2D NaN image)
%   ds - downsampling kernel (int, typically 32)
%   maxThreads - maximum number of cores to use (default 4)
%   initialParams - initial timing parameters to use (default defined 
%       below)
%   corrWin - window for performance over time [t1 t2]
% OUTPUTS:
%   perf - correlation map showing correlation between input HbT and
%       predicted HbT
%   IRF - estimated IRFs across the cortex
%   outParams - output parameters including t0, tauA, tauB and brain
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

ds_brain_mask = f_downsample(brain_mask,ds);
nanIdx = ~isnan(ds_brain_mask(:));
N = sum(nanIdx);

dim = size(dsHbT);
dsHbT = reshape(dsHbT,[],dim(3));
dsCa = reshape(dsCa,[],dim(3));

dsHbT = dsHbT(nanIdx,:)';
dsCa = dsCa(nanIdx,:)';

dsHbT = dsHbT./std(dsHbT,0,1);
dsCa = dsCa./std(dsCa,0,1);

%% create design matrices

n_IRF = fs*(win(2)-win(1))+1;

Ca_mat = zeros(dim(3),N,n_IRF);

T = size(dsCa,1);
l_irf = fs*range(win)+1;
idx_irf = win(1)*fs:win(2)*fs;

i1 = abs(min([idx_irf;zeros(1,numel(idx_irf))],[],1))+1;
i2 = T-max([idx_irf;zeros(1,numel(idx_irf))],[],1);
i3 = max([idx_irf;zeros(1,numel(idx_irf))],[],1)+1;
i4 = T-i1+1;

for v = 1:l_irf
    Ca_mat(i3(v):i4(v),:,v) = dsCa(i1(v):i2(v),:);
end

Ca_mat = permute(Ca_mat,[3 2 1]);
design_HbT = permute(dsHbT,[3 2 1]);

%% optimize initial parameters
if nargin > 8
    t0 = initialParam(1);
    tau1 = initialParam(2);
    tau2 = initialParam(3);
else
    t0 = 0.1;
    tau1 = 0.5;
    tau2 = 0.53;
end
A = ones(1,N);
B = -ones(1,N);

hrf1 = f_alpha_IRF(t0,tau1,tau2,A,0,fs,win);
hrf2 = f_alpha_IRF(t0,tau1,tau2,0,B,fs,win);

convPos = sum(Ca_mat.*hrf1,1);
convNeg = sum(Ca_mat.*hrf2,1);

LR = f_hemRegress(design_HbT,cat(4,convPos,convNeg),ones(1,N));

A = LR(:,:,1);
B = -LR(:,:,2);

%% run optimization using gradient descent
fun = @(params)f_hrf_cost_func(params(1),params(2),params(3),params(4:N+3),params(N+4:2*N+3),fs,win,Ca_mat,design_HbT);

options = optimset('MaxFunEvals',25000,'MaxIter',500,'Display','iter','Algorithm','active-set');
params = fmincon(fun,[t0, tau1, tau2, A, B],[],[],[],[],[0,0.01,0.01,-Inf*ones(1,numel(A)*2)],[win(2) win(2) win(2) Inf*ones(1,numel(A)*2)],[],options);

%% pixel-wise LR

hrf1 = f_alpha_IRF(params(1),params(2),params(3),1,0,fs,win);
hrf2 = f_alpha_IRF(params(1),params(2),params(3),0,-1,fs,win);

convPos = f_3Dconvolve(Ca./std(Ca,0,3),hrf1,win*fs,brain_mask);
convNeg = f_3Dconvolve(Ca./std(Ca,0,3),hrf2,win*fs,brain_mask);

LR = f_hemRegress(HbT./std(HbT,0,3),cat(4,convPos,convNeg),brain_mask);

outParams = struct;
outParams.t0 = params(1);
outParams.tauA = params(2);
outParams.tauB = params(3);
outParams.A = LR(:,:,1);
outParams.B = -LR(:,:,2);

%%
dim = size(Ca);

IRF = f_alpha_IRF(params(1),params(2),params(3),outParams.A(:)',outParams.B(:)',fs,win)';
IRF = reshape(IRF,dim(1),dim(2),[]);

%% calculate performance
pred_HbT = convPos.*LR(:,:,1)+convNeg.*LR(:,:,2);

perf = f_corr(HbT,pred_HbT,3);

if ~isempty(corrWin)
    corrGram = f_HemCorrGram(HbT,pred_HbT,corrWin);
end
%%

function J = f_hrf_cost_func(t0,tau1,tau2,A,B,sr,range,X_mat,y)
    hrf = f_alpha_IRF(t0,tau1,tau2,A,B,sr,range);
    
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