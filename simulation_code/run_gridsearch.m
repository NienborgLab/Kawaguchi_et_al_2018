function run_gridsearch
% simply run all the grid searches

models = {'ITB', 'ITBdt', 'Yates', 'Ralf', 'Attractor_Linear', 'Attractor_Sigmoid'};
lenm = length(models);
for m = 1:lenm
    disp([models{m} ' started...'])
    mdl = models{m};
    grid_search_like_cfnoise(mdl, 1);
end