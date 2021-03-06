function [pka_all, pka_hc, pka_lc, pkavar_all, pkavar_hc, pkavar_lc] = getPKA(ss, stm, ch, cf, nbin, pkmethod, repeat)
% compute psychophysical kernel amplitude (PKA) in a specified method
% INPUT:
% ss ... signed assigned stimulus (in disparity task, hdx x Dc)
% stm ... stimulus matrix (trial x frames)
% ch ... binary choice (0 or 1; trial x 1)
% cf ... confidence (trial x 1)
% nbin ... the number of time bins for PKA
% pkamethod ... method to compute PKA 
% (0; Nienborg & Cumming, 2009; 1: image classification, 2: logistic regression)
% repeat ... the number of repeat for resampling
%
% OUTPUT:
% pka_all ... averaged PKA
% pka_hc ... PKA in high confidence trials
% pka_lc ... PKA in low confidence trials
% pkavar_all ... variance of pka_all
% pkavar_hc ... variance of pka_hc
% pkavar_lc ... variance of pka_lc
%
% NOTE: when repeat > 0, resampling is performed and estimated SEM is 
% stored in the second row of each output pka
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% median split of confidence in each stimulus type
cf = binalize(cf, ss);

% compute PKA
[pka_all, pka_hc, pka_lc, pkavar_all, pkavar_hc, pkavar_lc] = ...
    computePKA(ss, stm, ch, cf, nbin, pkmethod);

% resampling for error bars
if nargin < 7; repeat = 0; end
if repeat > 0
    disp('resample')
    [err, err0, err1] = resamplePK(ss, stm, ch, cf, nbin, repeat, pkmethod);
    pka_all = [pka_all; err];
    pka_hc = [pka_hc; err1];
    pka_lc = [pka_lc; err0];
end

% subfunctions
function cf = binalize(cf, ss)
% median split of cf
u = unique(cf);
cf_orig = cf;
if length(u)==2
    cf(cf_orig==min(cf_orig)) = 0;
    cf(cf_orig==max(cf_orig)) = 1;
else
    unis = unique(ss);
    lens = length(unis);
    for s = 1:lens
        [cf0, cf1] = median_split(cf(ss==unis(s)));
        cf(ss==unis(s) & ismember(cf_orig, cf0)) = 0;
        cf(ss==unis(s) & ismember(cf_orig, cf1)) = 1;
    end
end

function [pka_all, pka_hc, pka_lc, pkavar_all, pkavar_hc, pkavar_lc] = ...
    computePKA(ss, stm, ch, cf, nbin, pkmethod)
switch pkmethod
    case 0 % Nienborg & Cumming, 2009
        [pka_all, pka_hc, pka_lc, pkavar_all, pkavar_hc, pkavar_lc] = PKA_hn(ss, stm, ch, cf, nbin);
    case 1 % image classification
        [pka_all, pka_hc, pka_lc, pkavar_all, pkavar_hc, pkavar_lc] = PKA_ic(ss, stm, ch, cf, nbin);
    case 2 % logistic regression
        [pka_all, pka_hc, pka_lc, pkavar_all, pkavar_hc, pkavar_lc] = PKA_logreg(stm, ch, cf, nbin);
end

function [cf0, cf1] = median_split(cf)
[~,idx] = sort(cf);
n = floor(length(cf)/2);
cf0 = cf(idx(1:n));
cf1 = cf(idx(n+1:end));

function stmbin = binstm(stm, nbin)
% bin the stimulus metrix along with time
nframe = size(stm, 2);
begin = 1;
frameperbin = floor(nframe/nbin);
stmbin = nan(size(stm, 1), nbin);
for a = 1:nbin
    stmbin(:, a) = mean(stm(:,begin:begin+frameperbin-1), 2);
    begin = begin + frameperbin;
end

function [pk, pkvar] = getKernel(stm, disval, ch)
% trial averaged stimulus distributions split by choice
nd = length(disval);
ntr = length(ch);
svmat = zeros(ntr, nd);
for r = 1:ntr
    for d = 1:nd
        svmat(r,d) = sum(stm(r,:)==disval(d));
    end
end
% compute PK for 0% stimulus
pk = mean(svmat(ch==1,:)) - mean(svmat(ch==0,:));
pkvar = var(svmat(ch==1, :), [], 1) + var(svmat(ch==0, :), [], 1);

function [pka_all, pka_hc, pka_lc, pkavar_all, pkavar_hc, pkavar_lc] = PKA_hn(ss, stm, ch, cf, nbin)
% only 0% signal trials
stm = stm(ss==0,:);
ch = ch(ss==0);
cf = cf(ss==0);
% time-averaged PK in each time bin
disval = unique(stm);
nd = length(disval);
% discretize stimuli, if too many unique values
if nd > 25
    [~,~,stm] = histcounts(stm, 11);
    stm = stm - mean(mean(stm));
end
disval = unique(stm);
nd = length(disval);
tkernel = nan(nd, nbin);
tkernel_h = nan(nd, nbin);
tkernel_l = nan(nd, nbin);
tkernelvar = nan(nd, nbin);
tkernelvar_h = nan(nd, nbin);
tkernelvar_l = nan(nd, nbin);
begin = 1;
nframe = size(stm, 2);
frameperbin = floor(nframe/nbin);
for a = 1:nbin
    [pk0, pk0var] = getKernel(stm(:,begin:begin+frameperbin-1), disval, ch);
    [pkh, pkhvar] = getKernel(stm(cf == 1, begin:begin+frameperbin-1), disval, ch(cf == 1));
    [pkl, pklvar] = getKernel(stm(cf == 0, begin:begin+frameperbin-1), disval, ch(cf == 0));
    tkernel(:,a) = pk0';
    tkernel_h(:,a) = pkh';
    tkernel_l(:,a) = pkl';
    tkernelvar(:,a) = pk0var';
    tkernelvar_h(:,a) = pkhvar';
    tkernelvar_l(:,a) = pklvar';
    begin = begin + frameperbin;
end
% PKA
% pk = getKernel(stm, disval, ch);
pk = mean(tkernel,2);
pkvar = mean(tkernelvar(:));
pka_all = nan(1,nbin);
pka_hc = nan(1,nbin);
pka_lc = nan(1,nbin);
pkavar_all = nan(1,nbin);
pkavar_hc = nan(1,nbin);
pkavar_lc = nan(1,nbin);
for a = 1:nbin
    pka_all(a) = dot(tkernel(:,a), pk);
    pka_hc(a) = dot(tkernel_h(:,a), pk);
    pka_lc(a) = dot(tkernel_l(:,a), pk);
    pkavar_all(a) = mean(tkernelvar(:,a),1)*pkvar;
    pkavar_hc(a) = mean(tkernelvar_h(:,a),1)*pkvar;
    pkavar_lc(a) = mean(tkernelvar_l(:,a),1)*pkvar;
end

function [pka_all, pka_hc, pka_lc, pkavar_all, pkavar_hc, pkavar_lc] = PKA_ic(ss, stm, ch, cf, nbin)
% only 0% signal trials
stm = stm(ss==0,:);
ch = ch(ss==0);
cf = cf(ss==0);
% image classification to compute PKA
pka_all = mean(stm(ch==1,:), 1) - mean(stm(ch==0,:), 1);
pka_hc = mean(stm(ch==1 & cf==1,:), 1) - mean(stm(ch==0 & cf==1,:), 1);
pka_lc = mean(stm(ch==1 & cf==0,:), 1) - mean(stm(ch==0 & cf==0,:), 1);
pka_all = binstm(pka_all, nbin);
pka_hc = binstm(pka_hc, nbin);
pka_lc = binstm(pka_lc, nbin);
pkavar_all = var(stm(ch==1,:), [], 1) + var(stm(ch==0,:), [], 1);
pkavar_hc = var(stm(ch==1 & cf==1,:), [], 1) + var(stm(ch==0 & cf==1,:), [], 1);
pkavar_lc = var(stm(ch==1 & cf==0,:), [], 1) + var(stm(ch==0 & cf==0,:), [], 1);

function [pka_all, pka_hc, pka_lc, pkavar_all, pkavar_hc, pkavar_lc] = PKA_logreg(stm, ch, cf, nbin)
% logistic regression (intercept included to capture a bias)
stm = zscore(binstm(stm, nbin));
pka_all = glmfit(stm, ch, ...
    'binomial', 'link', 'logit', 'constant', 'on');
pka_hc = glmfit(stm(cf == 1, :), ch(cf == 1), ...
    'binomial', 'link', 'logit', 'constant', 'on');
pka_lc = glmfit(stm(cf == 0, :), ch(cf == 0), ...
    'binomial', 'link', 'logit', 'constant', 'on');
pka_all = pka_all(2:end)'; pka_hc = pka_hc(2:end)'; pka_lc = pka_lc(2:end)';
pkavar_all = var(stm(ch==1,:), [], 1) + var(stm(ch==0,:), [], 1);
pkavar_hc = var(stm(ch==1 & cf==1,:), [], 1) + var(stm(ch==0 & cf==1,:), [], 1);
pkavar_lc = var(stm(ch==1 & cf==0,:), [], 1) + var(stm(ch==0 & cf==0,:), [], 1);

function [err, err0, err1] = resamplePK(ss, stm, ch, cf, nbin, repeat, pkmethod)
% resampling procedure
ampr = nan(repeat, nbin);
ampr0 = nan(repeat, nbin);
ampr1 = nan(repeat, nbin);
ntr = length(ch);
for r = 1:repeat
    disp(r)
    rtr = randi(ntr, ntr, 1);
    [ampr(r,:), ampr1(r,:), ampr0(r,:)] = computePKA(ss(rtr), stm(rtr,:), ch(rtr), cf(rtr), nbin, pkmethod);
end
err = std(ampr, [], 1);
err0 = std(ampr0, [], 1);
err1 = std(ampr1, [], 1);