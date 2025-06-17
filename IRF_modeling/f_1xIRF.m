function [perf,IRF,outParams,corrGram] = f_1xIRF(HbT,Ca,win,fs, ...
    brain_mask,ds,corrWin,maxThreads,initialParam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs deconvolution of HbT by applying a single double alpha function 
% convolution kernel to Ca. 
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
%   outParams - output parameters including t0, tauA, tauB, A, and B
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

Ca_mat = reshape(Ca_mat,[],n_IRF);
design_HbT = dsHbT(:);

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
A = 1;
B = -1;

hrf1 = f_alpha_IRF(t0,tau1,tau2,A,0,fs,win);
hrf2 = f_alpha_IRF(t0,tau1,tau2,0,B,fs,win);

convPos = sum(Ca_mat.*hrf1',2);
convNeg = sum(Ca_mat.*hrf2',2);

LR = [convPos,convNeg] \ design_HbT;
A = LR(1);
B = -LR(2);

%% run optimization using gradient descent
fun = @(params)f_hrf_cost_func(params(1),params(2),params(3),params(4),params(5),fs,win,Ca_mat,design_HbT);

options = optimset('MaxFunEvals',25000,'MaxIter',500,'Display','iter','Algorithm','active-set','FunValCheck','on');
params = fmincon(fun,[t0, tau1, tau2, A, B],[],[],[],[],[0,0.01,0.01,-Inf,-Inf],[win(2),win(2),win(2),Inf,Inf],[],options);

%%

IRF = f_alpha_IRF(params(1),params(2),params(3),params(4),params(5),fs,win);
conv = f_3Dconvolve(Ca./std(Ca,0,3),IRF,win*fs,ones(size(brain_mask)));
perf = f_corr(HbT,conv,3);

if ~isempty(corrWin)
    corrGram = f_HemCorrGram(HbT,conv,corrWin);
end
%%

outParams = struct;
outParams.t0 = params(1);
outParams.tau1 = params(2);
outParams.tau2 = params(3);
outParams.A = params(4);
outParams.B = params(5);

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = f_hrf_cost_func(t0,tau1,tau2,A,B,sr,range,X_mat,y)
    hrf = f_alpha_IRF(t0,tau1,tau2,A,B,sr,range);
    
    conv_result = sum(X_mat.*hrf',2);
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
