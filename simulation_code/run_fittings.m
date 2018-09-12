function run_fittings(animal)
% simply run all the fittings

models = {'ITB', 'ITBdt', 'Yates', 'Ralf', 'Attractor_Linear', 'Attractor_Sigmoid'};
nstart = [1, 1, 1, 1, 1, 1];
lenm = length(models);
parfor m = 1:lenm
    mdl = models{m};
    fitmodel2pka(mdl, animal, nstart(m), 1);
end
