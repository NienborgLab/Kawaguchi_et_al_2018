function results = grid_search_like(model, saveoption)
%% 
% perform a grid search on specified model parameters 
% This function is to explore how PKA changes based on each model's
% parameters
%
% INPUT: 
% model ... 'ITB', 'ITBdt', 'Yates', 'Ralf', 'Attractor_Linear', 'Attractor_Sigmoid'
% saveoption ... 0, no save 1, save
%
% EXAMPLE: results = grid_search_like('ITB', 0);
%

% path
mypath = pathfinder;

%%
% dealing with inputs
if nargin < 2; saveoption = 0; end
switch model
    case {'ITB',  'ITBdt'}
        % decision bound
        params{1} = [0:0.5:20]*0.2591;
        % combination
        idxmat = combvec(1:length(params{1}));
    case 'Race'
        % decision bound
        params{1} = [0:1:80]*16.5831;
        % combination
        idxmat = combvec(1:length(params{1}));
    case 'Yates'
        % tmax
        params{1} = 100:100:500;
        % tau
        params{2} = 10:20:100;
        % offset_gain
        params{3} = 0:0.3:1.5;
        % kernel gain ratio
        params{4} = 0:0.25:3;
        % combination
        idxmat = combvec(1:length(params{1}), ...
            1:length(params{2}), 1:length(params{3}), ...
            1:length(params{4}));
    case {'Attractor_Linear', 'Attractor_Sigmoid'}
        % time
        params{1} = 20:20:1000;
        % acceleration
        params{2} = 0.01:0.01:0.4;
        % combination
        idxmat = combvec(1:length(params{1}), ...
            1:length(params{2}));
    case 'Ralf'
        % time
        params{1} = 12:3:99;
        % combination
        idxmat = combvec(1:length(params{1}));        
end

results.model = model;
results.params = params;

%%
% grid search
niter = size(idxmat, 2);
npara = length(params);
disp([model ' model: '])
disp([num2str(npara) ' parameters: ' num2str(niter) ' iterations'])
iterparams = zeros(niter, npara);
for n = 1:niter
    for k = 1:npara
        iterparams(n, k) = ...
            params{k}(idxmat(k,n));
    end
    results.iteration(n).params = iterparams(n,:);
end

temp = cell(niter, 1);
disp('started grid search...')
switch model
    case 'ITB'
        parfor n = 1:niter            
            temp{n} = SDTfast('ntr',50000, 'db', iterparams(n, 1));
        end
    case 'ITBdt'
        parfor n = 1:niter
            temp{n} = SDTfast('ntr',50000, 'db', iterparams(n, 1), 'dt');
        end
    case 'Race'
        parfor n = 1:niter
            temp{n} = SDT_PKA('ntr',100000, 'race', 'db', iterparams(n, 1), ...
                'conftype','Bayes','pkmethod',0);
        end
    case 'Yates'
        parfor n = 1:niter
            p = iterparams(n,:);
            out = Yates_simulation('ntr',5000,'tmax',p(1),...
                'tau', p(2), 'offset_gain', p(3), ...
                'kernelgain_s', 0.05, 'kernelgain_c', 0.05*p(4));
            temp{n} = struct('psth', out.psth, 'pka', out.pka, ...
                'pka_highconf', out.pka_highconf, 'pka_lowconf', ...
                out.pka_lowconf);
        end
    case 'Attractor_Linear'
        parfor n = 1:niter
            p = iterparams(n,:);
            temp{n} = EvidenceAccumulation(p(2), 1, ...
                p(1), 50000, 'normal', 'linear', 0, 0);
        end
    case 'Attractor_Sigmoid'
        parfor n = 1:niter
            p = iterparams(n,:);
            temp{n} = EvidenceAccumulation(p(2), 1, ...
                p(1), 50000, 'normal', 'sigmoid', 0, 0);
        end
    case 'Ralf'
        % load output structure
        try
            load([mypath '/data/E.mat'])
        catch
            error('No output struct  E is found. To obtain the output struct E, follow the instruction in https://github.com/katsu1110/sampling_decision')
        end
        for n = 1:niter
            p = iterparams(n,:);
            [pkas, pka0] = Kernel_Compute(trcut(E, p(1)));
            temp{n} = [pka0; flipud(pkas)];
        end        
end

% simulation results storage
disp('finished a grid search. Now storing data into struct...')
switch model
    case {'ITB', 'ITBdt'}
        for n = 1:niter
           results.iteration(n).pka =  [temp{n}.amplitude; ...
                    temp{n}.amplitude_highconf; temp{n}.amplitude_lowconf];      
           results.iteration(n).noise = [temp{n}.noisestm, temp{n}.noiseidv];
        end
    case 'Race'
        for n = 1:niter
           results.iteration(n).pka =  [temp{n}.pka; ...
                    temp{n}.pka_highconf; temp{n}.pka_lowconf];      
           results.iteration(n).noise = [temp{n}.noisestm, temp{n}.noiseidv];
        end
    case 'Yates'
        for n = 1:niter
            results.iteration(n).psth = temp{n}.psth;
            results.iteration(n).pka = [temp{n}.pka; temp{n}.pka_highconf; temp{n}.pka_lowconf];
        end
    case {'Attractor_Linear', 'Attractor_Sigmoid', 'Ralf'}
        for n = 1:niter
            results.iteration(n).pka = temp{n};
        end
end
disp('done')

% autosave
if saveoption==1
    save([mypath '/data/struct_storage/grid_' model '.mat'], 'results')
    disp('saved!')
end