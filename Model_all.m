function Model_all(option)
% reproduce main results in Katsuhisa et al., 2018
%
% INPUT:
% option ... 0, start over from a 'grid search'
%            1, use already stored data 
%
% EXAMPLE: Model_all(1)
%

% path
mypath = cd;
addpath(genpath(mypath))

% inputs
if nargin < 1; option = 0; end

%% 
% signatures of confidence (prediction of SDT) (Figure 5b)
figure;
[cf, acc, stm] = confidence_signature(100000);
plot_confidence_signature(cf, acc, stm);
set(gcf, 'Name', 'signature of confidence by SDT', 'NumberTitle', 'off')
try
    savefig([mypath '/data/figure_storage/sdt_confsig.fig'])
catch
    disp('saving figure is skipped due to no data/figure_storage folder')
end
disp('saved!')
disp('--------------------------------------------------------')

%%
% cost on having a decision bound on performance (Figure 9)
figure;
db_cost
set(gcf, 'Name', 'cost of having DB', 'NumberTitle', 'off')
try
    savefig([mypath '/data/figure_storage/dbcost.fig'])
catch
    disp('saving figure is skipped due to no data/figure_storage folder')
end
disp('saved!')
disp('--------------------------------------------------------')

%% 
% SDT-based model prediction (Hangya et al., 2016) to the third confidence
% signature (used in Figure 6d)
try
    sdt_prediction_animal
catch
    disp('saving figure is skipped due to no data/figure_storage folder')
end
disp('--------------------------------------------------------')

%%
% 'grid search' (Figure 1, 8)
models = {'ITB', 'ITBdt', 'Yates', 'Attractor_Linear', 'Attractor_Sigmoid', 'Ralf'};
lenm = length(models);
if option == 0
    disp('start a grid search ... ')
    for m = 1:lenm
        grid_search_like(models{m}, 1);
    end
end
for m = 1:lenm
    figure;
    try
        plot_gridsearch(models{m})
    catch
        disp('saving figure is skipped due to no data/figure_storage folder')
    end
end
disp('--------------------------------------------------------')

%% 
% model fitting (Figure 7; models)
if option == 0
    disp('start a model fitting ...')
    disp('Animal A...')
    run_fittings('kiwi')
    disp('Animal B...')
    run_fittings('mango')
end
figure;
try
    fitted_modelcompare
catch
    disp('saving figure is skipped due to no data/figure_storage folder')
end
disp('--------------------------------------------------------')