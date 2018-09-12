function fitted_modelcompare
% perform model comparison based on each model's fit to pka data

% path
mypath = pathfinder;

% yellow and green
yellow = [0.9576    0.7285    0.2285];
green = [0.1059    0.4706    0.2157];

% models to compare with
models = {'ITB','ITBdt','Yates','Ralf'};
modelnames = {'ITB','ITBdt','Yates','Ralf'};
lenm = length(models);
close all;
animals = {'mango', 'kiwi'};
lena = length(animals);
mse = nan(lena, lenm);
aic = nan(lena, lenm);
nobs = 8;
for m = 1:lenm
    disp([models{m} ' is on process...'])    
    for a = 1:lena        
        % load fitting results
        load([mypath '/data/struct_storage/pkafit_' models{m} '_' animals{a} '.mat']);
        
        % animals' pka
        switch animals{a}
            case 'mango'
                variance = [0.17858 0.16075 0.15057 0.14073, ...
                    0.15405 0.15508 0.15226 0.15274];
            case 'kiwi'
                variance = [0.033055 0.032098 0.028395 0.02767, ...
                    0.038769 0.039468 0.040962 0.040099];
        end

        % recompute mse, nmse, AIC & BIC properly
        beta = glmfit(results.pred', results.data', 'normal', 'link', 'identity', 'constant', 'off');
        results.mse = sum((beta(1)*results.pred - results.data).^2)/nobs;
        results.nmse = sum(((beta(1)*results.pred - results.data).^2)./variance)/nobs;
        L = -(nobs/2)*log(2*pi*mean(variance)) - nobs*results.mse/(2*mean(variance));
        results.loglikelihood = L;
        [results.aic, results.bic] = aicbic(L, length(results.params), nobs);
        
        % store values
        save([mypath '/data/struct_storage/pkafit_' models{m} '_' animals{a} '.mat'], 'results');
        mse(a,m) = results.nmse;
        aic(a,m) = results.aic;
    end
end
lenmet = 2;
for a = 1:lena
    
    % MSE
    subplot(lenmet,lena,a)
    bar(1:lenm, mse(a,:),'facecolor','w','edgecolor','k')
    ylabel({'normalized','MSE'})
    set(gca,'XTick',1:lenm,'XTickLabel',{'','','','','',''})
    xtickangle(45)
    set(gca,'box','off'); set(gca,'TickDir','out');
    title(animals{a})
    
    % AIC
    subplot(lenmet,lena,a+lena)
    mina = min(aic(a,:));
    bar(1:lenm, aic(a,:) - mina,'facecolor','w','edgecolor','k')
    ylabel('\Delta AIC')
    set(gca,'XTick',1:lenm,'XTickLabel',modelnames)
    xtickangle(45)
    set(gca,'box','off'); set(gca,'TickDir','out');
        
%     % exact values
%     for m = 1:lenm
%         % add values
%         subplot(lenmet,lena,a)
%         text(m-0.55, mse(a,m), num2str(mse(a,m)))
%         subplot(lenmet,lena,a+lena)
%         text(m-0.55, aic(a,m) - mina, num2str(aic(a,m)-mina))
%     end    
end
set(gcf, 'Name', 'modelcompare_metric','NumberTitle','off')
savefig([mypath '/data/figure_storage/modelcompare_metric.fig']);

% best fit in each model
h = figure;
for a = 1:lena
    for m = 1:lenm
        subplot(lena, lenm, m+(a-1)*lenm)
        % load fitting results
        load([mypath '/data/struct_storage/pkafit_' models{m} '_' animals{a} '.mat']);
        % rescaling
        beta = glmfit(results.pred', results.data', 'normal', 'link', 'identity', 'constant', 'off');
        % model prediction
        p = plot(1:4, beta*results.pred(5:8), '-', 'color', green, 'linewidth', 2);
        p.Color(4) = 0.4;
        hold on;
        p = plot(1:4, beta*results.pred(1:4), '-', 'color', yellow, 'linewidth', 2);
        p.Color(4) = 0.4;
        hold on;
        % data
        scatter(1:4, results.data(5:8), 20, 'filled', 'markerfacecolor', green)
        hold on;
        scatter(1:4, results.data(1:4), 20, 'filled', 'markerfacecolor', yellow)
        set(gca,'XTick',1:lenm,'XTickLabel',modelnames)
        xlim([0.75 4.25])
        if a==1
            title(modelnames{m})
        end
        set(gca, 'XTick', [1 4], 'XTickLabel', [0 1500])
        if m==1
            ylabel({'normalized','PKA',animals{a}})
            if a==2
                xlabel('time from stimulus onset (ms)')
            end
        end
%         % parameters
%         xx = get(gca, 'XLim');
%         yy = get(gca, 'YLim');
%         for n = 1:length(results.params)
%             % normalized decision bound
%             if m < 3 && n==2
%                 results.params(n) = results.params(n)/0.2591;
%             end
%             text(xx(1)+0.7*(xx(2)-xx(1)), yy(1)+(1 - 0.1*n)*(yy(2)-yy(1)), ...
%                 num2str(results.params(n)))
%         end
        set(gca,'box','off'); set(gca,'TickDir','out');
    end
end
set(h, 'Name', 'modelcompare_bestfits', 'NumberTitle', 'off')
savefig([mypath '/data/figure_storage/modelcompare_bestfits.fig']);
disp('figure saved')