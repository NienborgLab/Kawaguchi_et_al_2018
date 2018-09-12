function sdtpred = sdt_prediction_animal
% SDT prediction for the third confsig

% path
mypath = pathfinder;

% data from animals
animals = {'Mango', 'Kiwi'};
noise = [45.5956, 43.9113; 23.92, 19.4023];
for v = 1:2    
    for a = 1:2
        % animals' data
        sdtpred.animal(a).avrew(v).noise = noise(a, v);
        switch a
            case 1 % Mango
                sdtpred.animal(a).avrew(v).dc = 100*[0    0.0300    0.0600    0.1250    0.2500    0.5000    1.0000];                
            case 2 % Kiwi
                sdtpred.animal(a).avrew(v).dc = 100*[0    0.0310    0.0625    0.1250    0.2500    0.5000];
        end
        sdtpred.animal(a).avrew(v).name = animals{a};   
        if v==1 && a==1
            cor = [0.4927    0.4738    0.4752    0.5116    0.5687    0.6666    0.7405];
            err = [0.4925    0.4427    0.4496    0.4729];
        elseif v==1 && a==2
            cor = [0.3042 0.3423 0.3563 0.3914 0.5031 0.7397];
            err = [0.2959 0.2906 0.2809];
        elseif v==2 && a==1
            cor = [0.5716 0.5921 0.5986 0.6229 0.7002 0.8369 0.9012];
            err = [0.6264 0.5734 0.5872 0.5996];
        elseif v==2 && a==2
            cor = [0.4130 0.4502 0.4779 0.5277 0.6940 0.9972];
            err = [0.4516 0.4270 0.4283];
        end
        sdtpred.animal(a).avrew(v).data.correct_conf = cor;
        sdtpred.animal(a).avrew(v).data.error_conf = err;
        
        % SDT prediction
        [cf, acc, stm] = confidence_signature(100000, sdtpred.animal(a).avrew(v).noise, ...
            'sdt', 'uniform', sdtpred.animal(a).avrew(v).dc);
        sdtpred.animal(a).avrew(v).sdt_pred = plot_confidence_signature(cf, acc, stm, 0);
    end
end

% autosave
save([mypath '/data/struct_storage/sdtpred.mat'], 'sdtpred')
disp('saved!')