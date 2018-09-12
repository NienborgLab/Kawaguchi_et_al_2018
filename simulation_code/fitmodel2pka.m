function results = fitmodel2pka(model, animal, nstart, saveoption)
%% 
% fit a specified model to the animals' data (PKA). Confidence noise is
% added on top of each model's free parameters to represent imperfect
% relationship between confidence and proportion correct.
%
% INPUT: 
% model ... 'ITB', 'ITBdt', 'Yates', 'Ralf', 'Attractor_Linear', 'Attractor_Sigmoid'
% animal ... 'kiwi' or 'mango'
% nstart ... the number of starting points 
% saveoption ... 0, no save 1, save
%
% EXAMPLE: results = fitmodel2pka('ITB', 'mango', 0);
%

%%
% path
mypath = pathfinder;

% dealing with inputs
if nargin < 1; animal = 'ITB'; end
if nargin < 2; animal = 'kiwi'; end
if nargin < 3; nstart = 1; end
if nargin < 4; saveoption = 0; end

% initial parameters
switch model
    case {'ITB',  'ITBdt'}
        % decision bound
        p0 = [0, 2.4];
        lbd = [0, 1*0.2591];
        hbd = [3, 15*0.2591];
    case 'Race'
        % decision bound
        p0 = [0, 27*16.5831];
        lbd = [0, 5*16.5831];
        hbd = [3, 50*16.5831];
    case 'Yates'
        % tmax, tau, offset_gain, kernel gain ratio
        p0 = [2.5, 420, 24, 0.3, 1.5];
        lbd = [0, 100, 10, 0, 0];
        hbd = [3, 500, 100, 1.5, 3];
    case 'Attractor_Linear'
        % trial duration, acceleration
        p0 = [0, 40, 0.04];
        lbd = [0, 20, 0.01];
        hbd = [3, 1000, 0.4];
    case 'Attractor_Sigmoid'
        % trial duration, acceleration
        p0 = [2.9, 100, 0.03];
        lbd = [0, 20, 0.01];
        hbd = [3, 1000, 0.4];    
    case 'Ralf'
        % trial duration
        p0 = [1.2, 30];
        lbd = [0 12];
        hbd = [3 99];
        % load output struct
        try
            load([mypath '/data/E.mat'])
        catch
            error('No output struct  E is found. To obtain the output struct E, follow the instruction in https://github.com/katsu1110/sampling_decision')
        end
end

% animals' pka
switch animal
    case 'mango'
        nobs = 8;
        variance = [0.17858 0.16075 0.15057 0.14073, ...
                    0.15405 0.15508 0.15226 0.15274];
        data = [1.2965, 0.9326, 0.7356, 0.3342, ...
            0.7035, 0.7095, 0.6630, 0.4813];
    case 'kiwi'
        nobs = 8;
        variance = [0.033055 0.032098 0.028395 0.02767, ...
                    0.038769 0.039468 0.040962 0.040099];
        data = [1.3224, 1.1016, 0.6477, 0.4618, ...
            0.6776, 0.6542, 1.0576, 0.7000];
end

% cost function & parameter estimate
switch model
    case 'Ralf'
        c = @(p)cost_Ralf(p, E, data);
    otherwise
        c = @(p)cost(p, model, data);
end
options = optimset('MaxFunEvals',10000,'maxiter',10000);
lenp = length(p0);
params = zeros(nstart, lenp);
for n = 1:nstart
    disp(['fitting started: trial ' num2str(n) '/' num2str(nstart)])
    p_init = p0;
    for p = 1:lenp
        p_init(p) = normrnd(p_init(p), 1);
        if p_init(p) < lbd(p)
            p_init(p) = lbd(p);
        elseif p_init(p) > hbd(p)
            p_init(p) = hbd(p);
        end
    end
    params(n, :) = fminsearchbnd(c, p_init, lbd, hbd, options);
end
results.params_all = params;
params = mean(params, 1);

% into struct
results.model = model;
results.animal = animal;
results.data = data;
% estimated parameters and prediction
results.params = params;
switch model
    case 'Ralf'
        tcand = 12:3:99;
        [~, idx] = min(abs(tcand - params(2)));
        [pkas, pka0] = Kernel_Compute(trcut(E, tcand(idx)), 'nbin', 4, 'cfnoise', params(1));
        results.pred = [pkas(2,:), pkas(1,:)]/max(pka0);
    otherwise
        results.pred = model_prediction(model, params);
end
% rescaling
beta = glmfit(results.pred', results.data', 'normal', 'link', 'identity', 'constant', 'off');
% mse
results.mse = sum((results.data - beta(1)*results.pred).^2)/nobs;
% chi
results.chi = sum(((results.data - beta(1)*results.pred).^2)./variance);
% log likelihood
results.loglikelihood = -(nobs/2)*log(2*pi*mean(variance)) - nobs*results.mse/(2*mean(variance));
% correlation
[rr, pp] = corrcoef(results.pred', results.data');
results.pearson_r = rr(1,2);
results.pearson_p = pp(1,2);
% AIC & BIC
[results.aic, results.bic] = aicbic(results.loglikelihood, length(params), nobs);

% autosave
if saveoption==1
    save([mypath '/data/struct_storage/pkafit_' model '_' animal '.mat'], 'results')
    disp('saved!')
end

%%
% subfunction
function pred = model_prediction(model, params)
switch model
    case 'ITB'
        temp = SDTfast('ntr',50000, 'db', params(2), 'cfnoise', params(1));
        pred = [temp.amplitude_highconf, temp.amplitude_lowconf]...
            /max(temp.amplitude);
    case 'ITBdt'
        temp = SDTfast('ntr',50000, 'db', params(2), 'dt', 'cfnoise', params(1));
        pred = [temp.amplitude_highconf, temp.amplitude_lowconf]...
            /max(temp.amplitude);
    case 'Race'
        temp = SDT_PKA('ntr',100000, 'race', 'db', params(2),...
            'cfnoise', params(1),'conftype','Bayes','pkmethod',0);
        pred = [temp.pka_highconf, pka.amplitude_lowconf]...
            /max(pka.amplitude);
    case 'Yates'
        temp = Yates_simulation('ntr',50000,'tmax',params(2),...
            'tau', params(3), 'offset_gain', params(4), ...
            'kernelgain_s', 0.05, 'kernelgain_c', 0.05*params(5), ...
            'cfnoise', params(1));
        pred = [temp.pka_highconf, temp.pka_lowconf]/max(temp.pka);   
        % remove biologically implausible parameter sets
        adap = temp.psth.stm_pref_adaptation(4);
        cp = pi*(mean(temp.psth.cp(1,2:end)) - 0.5)/sqrt(2);
        if adap < 0 || cp < 0
            pred = zeros(1,8);
        end
    case 'Attractor_Linear'
        temp = EvidenceAccumulation(params(3), 1, ...
            round(params(2)), 100000, 'normal', 'linear', params(1), 0);
        temp = tcbin(temp, 4);
        pred = [temp(3,:), temp(2,:)]/max(temp(1,:));
    case 'Attractor_Sigmoid'
        temp = EvidenceAccumulation(params(3), 1, ...
            round(params(2)), 100000, 'normal', 'sigmoid', params(1), 0);
        temp = tcbin(temp, 4);
        pred = [temp(3,:), temp(2,:)]./max(temp(1,:));
end

% cost function
function c = cost(p, model, data)
% model prediction
pred = model_prediction(model, p);

% rescaling
beta = glmfit(pred', data', 'normal', 'link', 'identity', 'constant', 'off');

% mse
c = sum((data - beta(1)*pred).^2)/8;

function c = cost_Ralf(p, E, data)
% model prediction
tcand = 12:3:99;
[~, idx] = min(abs(tcand - p(2)));
[pkas, pka0] = Kernel_Compute(trcut(E, tcand(idx)), 'nbin', 4, 'cfnoise', p(1));
pred = [pkas(2,:), pkas(1,:)]/max(pka0);

% rescaling
beta = glmfit(pred', data', 'normal', 'link', 'identity', 'constant', 'off');

% mse
c = sum((data - beta(1)*pred).^2)/8;