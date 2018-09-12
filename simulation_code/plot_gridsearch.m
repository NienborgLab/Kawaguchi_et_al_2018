function plot_gridsearch(model)
%% 
% Visualize results from 'grid_search_like.m' simulations
%
% INPUT: model ... 'ITB', 'ITBdt', 'Yates', 'Ralf', 'Attractor_Linear',
% 'Attractor_Sigmoid'
%
% EXAMPLE: plot_gridsearch('ITB')
%

% path
mypath = pathfinder;

% load struct
disp('loading data...')
load([mypath '/data/struct_storage/grid_' model '.mat'])
disp('loaded!')

% yellow and green
yellow = [0.9576    0.7285    0.2285];
green = [0.1059    0.4706    0.2157];

% close all;
boi = 4;
h = figure;
niter = length(results.iteration);
npara = length(results.params);
x = cell(1, npara);
d = zeros(1, niter);
hv = zeros(1, niter);
lv = zeros(1, niter);
switch model
    case {'ITB', 'ITBdt', 'Race'}
        for i = 1:niter
            x{1}(i) = results.iteration(i).params(1)...
                /results.iteration(i).noise(2);
            [d(i), hv(i), lv(i)] = normPKAd(results.iteration(i).pka, 0);
        end
        xlab = 'decision bound';
    case 'Yates'
        r = zeros(1, niter);
        l = zeros(1, niter);
        for i = 1:niter
            x{1}(i) = results.iteration(i).psth.stm_pref_adaptation(boi);
            x{2}(i) = pi*(results.iteration(i).psth.cp(boi) - 0.5)/sqrt(2);
            x{3}(i) = pi*(results.iteration(i).psth.cp(1) - 0.5)/sqrt(2);
            x{4}(i) = pi*(mean(results.iteration(i).psth.cp) - 0.5)/sqrt(2);
            d(i) = normPKAd(results.iteration(i).pka, 0);
            rr = corrcoef(1:size(results.iteration(i).pka,2), ...
                results.iteration(i).pka(1,:)/max(results.iteration(i).pka(1,:)));
            r(i) = rr(1,2);
            beta = glmfit(1:size(results.iteration(i).pka,2), ...
                results.iteration(i).pka(1,:)/max(results.iteration(i).pka(1,:)),...
                'normal','link','identity');
            l(i) = beta(2);
        end
        
        out = x{1} < -1 | isnan(d) | abs(x{2}) >= 1 | abs(x{3}) >= 1 | abs(x{4}) > 1;
        x{1}(out) = [];
        x{2}(out) = [];
        x{3}(out) = [];
        x{4}(out) = [];
        d(out) = [];
        r(out) = [];
        l(out) = [];
        [d, sidx] = sort(d);
        x{1} = x{1}(sidx);
        x{2} = x{2}(sidx);
        x{3} = x{3}(sidx);
        x{4} = x{4}(sidx);
        r = r(sidx);
        l = l(sidx);
        xlab = 'adaptation';
        ylab = 'CC';
    case {'Attractor_Linear', 'Attractor_Sigmoid'}
        for i = 1:niter
            x{1}(i) = results.iteration(i).params(1);
            x{2}(i) = results.iteration(i).params(2);
            norm = mean(results.iteration(i).pka(1, 1:2));
            hv(i) = mean(results.iteration(i).pka(3,end-1:end))/norm;
            lv(i) = mean(results.iteration(i).pka(2,end-1:end))/norm;
            d(i) = hv(i) - lv(i);
        end
        xlab = 'acceleration';
        ylab = 'time (a.u.)';
    case 'Ralf'
        for i = 1:niter
            x{1}(i) = results.iteration(i).params(1);
            results.iteration(i).pka(:,end) = results.iteration(i).pka(:,end-1);
            [d(i), hv(i), lv(i)] = normPKAd(results.iteration(i).pka, 0);
        end
        xlab = 'time (a.u.)';
end

switch model
    case 'Yates'    
        x{1}(r >= 0 | l >= 0) = [];
        x{2}(r >= 0 | l >= 0) = [];
        x{3}(r >= 0 | l >= 0) = [];
        x{4}(r >= 0 | l >= 0) = [];
        d(r >= 0 | l >= 0) = [];
        % insets       
        results.iteration(out) = [];
        results.iteration = results.iteration(sidx);
        results.iteration(r >= 0 | l >= 0) = [];
        t = [10 1224 753]; select = t;
        for s = 1:length(t)
            subplot(2,9, [9+s*2 10+s*2])
            pka = results.iteration(t(s)).pka(1,:);
            pkal = results.iteration(t(s)).pka(3,:);
            pkah = results.iteration(t(s)).pka(2,:);
            norm = max(pka(1,:));
            plot([1 size(pka,2)], [0 0], ':', 'color', 'k')
            hold on;
            plot(1:size(pkal,2), pkal/norm, '-', 'color', green)
            hold on;
            plot(1:size(pkah,2), pkah/norm, '-', 'color', yellow)
            xlim([1 size(pka,2)])
            ylim([-0.5 2])
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        end
        for k = 1:3
            switch k
                case 1
                    y = x{2};
                    nb = 1:4;
                case 2
                    y = x{3};
                    nb = 1;
                case 3
                    y = x{4};
                    nb = 4;
            end            
            subplot(2,9,[(k-1)*3+1 (k-1)*3+2 (k-1)*3+3])            
            scatter(x{1},y,[],d, 'filled', 'markerfacealpha', 0.2, 'markeredgealpha', 0.4)            
            hold on;                    
            for s = 1:length(select)
                plot(x{1}(select(s)),y(select(s)),'o','color',[0.8 0.6 1],'linewidth',1)
                hold on;
            end
            xx = [-1 1];
            yy = [-1 1];
            plot(xx, [0 0], ':k')
            hold on;
            plot([0 0], yy, ':k')
            hold on;
            plot(0.62*[1 1], yy, ':', 'color',0.5*[1 1 1],'linewidth',1.5)
            axis([xx yy])         
            c = colorbar;
            caxis([-1 1])
            c.Label.String = ['\Delta PKA (high - low conf) at time ' num2str(nb)];
            c.Label.Rotation = -90;
            xlabel(xlab)
            ylabel(ylab)
        end              
        
    case {'Attractor_Linear', 'Attractor_Sigmoid'}     
            % representitive params
            if strcmp(model, 'Attractor_Linear')
                acc = 0.1;
                trdur = 100;
                ext = [100 100 100; 0.01, 0.04, 0.1];
                xrange = [100, 0.1];
            elseif strcmp(model, 'Attractor_Sigmoid')
                acc = 0.06;
                trdur = 1000;
                ext = [20, 380, 1000; 0.06, 0.06, 0.06];
                xrange = [1000, 0.4];
            end        
            
            % vary trial duration
           subplot(2,4,1)
           idx = find(abs(x{2} - acc) < 0.005);
           scatter(x{1}(idx), hv(idx), 50, 'filled',...
               'markerfacecolor', yellow, 'markerfacealpha', 0.4)
           hold on;
           scatter(x{1}(idx), lv(idx), 50, 'filled',...
               'markerfacecolor', green, 'markerfacealpha', 0.4)
            xlim([min(x{1}) xrange(1)])
            xlabel('trial duration (a.u.)')
            ylabel('PKA_{t_{last}}');
            
           % vary acceleration
           subplot(2,4,2)
           idx = find(abs(x{1} - trdur) < 10);
           scatter(x{2}(idx), hv(idx), 50, 'filled',...
               'markerfacecolor', yellow, 'markerfacealpha', 0.4)
           hold on;
           scatter(x{2}(idx), lv(idx), 50, 'filled',...
               'markerfacecolor', green, 'markerfacealpha', 0.4)
           xlabel('acceleration')
           xlim([min(x{2}) xrange(2)])
            
           % acceleration vs trial duration 
           subplot(2,4,[3 4])
           Z = reshape(d, [length(unique(x{1})), ...
               length(unique(x{2}))]); 
            imagesc(unique(x{1}), ...
               unique(x{2}), Z);
            colormap(jet)
            caxis([-1 1.4])
            colorbar('eastoutside')
            set(gca,'XTick', unique(unique(x{1})))
            set(gca,'YTick', unique(x{2}))
            xlabel(ylab);
            ylabel(xlab);
            
            % inset
            for s = 1:3
                subplot(2,4,4+s)
                c1 = find(abs(x{1} - ext(1,s)) < 10);
                c2 = find(abs(x{2} - ext(2,s)) < 0.005);
                pka = results.iteration(intersect(c1, c2)).pka;
                [pka1] = tcbin(pka(1,:), ext(1,s)/2);
                [pka2] = tcbin(pka(2,:), ext(1,s)/2);
                [pka3] = tcbin(pka(3,:), ext(1,s)/2);
                norm = max(pka1);
                plot([1 size(pka1,2)], [0 0], ':', 'color', 'k')
                hold on;
                plot(1:size(pka1,2), pka2/norm, '-', 'color', green)
                hold on;
                plot(1:size(pka1,2), pka3/norm, '-', 'color', yellow)
                title({['trial duration: ' num2str(ext(1, s))], ['acceleration', num2str(ext(2, s))]})
                xlim([1 size(pka1,2)])
                ylim([-0.5 2])
            end
    case {'ITB','ITBdt','Race', 'Ralf'}
        subplot(2,5,[1 2 3 4])
        scatter(x{1}(:), hv, 60,  'markerfacecolor', yellow, ...
            'markerfacealpha', 0.4, 'markeredgecolor', 'w', 'markeredgealpha', 0.8)   
        hold on;
        scatter(x{1}(:), lv, 60, 'markerfacecolor', green, ...
            'markerfacealpha', 0.4, 'markeredgecolor', 'w', 'markeredgealpha', 0.8)      
        legend('high','low','location','northeast')
        legend('boxoff')
        legend('autoupdate','off')
        xx = get(gca, 'XLim');
        hold on;
        plot(xx, [0 0], ':k')
        xlabel(gca, xlab)
        ylabel('PKA_{tlast}')
        xlim([min(x{1}) max(x{1})])
        set(gca, 'XTick', [x{1}(1) x{1}(end)], 'XTickLabel', [x{1}(1) round(x{1}(end))])
        
        % insets
        db = arrayfun(@(x) x.params(1), results.iteration);
        unidb = unique(db);        
        noise = 0.2591;
        switch model
            case 'ITB'
                k = [unidb(9), unidb(15),  unidb(21), unidb(24), unidb(29)];
            case 'ITBdt'
                k = [unidb(9), unidb(15),  unidb(20), unidb(24), unidb(29)];
            case {'Race', 'Racedt'}
                k = [unidb(7), unidb(21),  unidb(32), unidb(40), unidb(48)];
%                 k = [unidb(7), unidb(16),  unidb(31), unidb(40), unidb(48)];
                noise = 16.5831;
            case 'Ralf'
                k = [unidb(5), unidb(9),  unidb(16), unidb(24), unidb(30)];
                noise = 1;
        end
        for s = 1:length(k)
            subplot(2,5, [5 + s])
            pka = results.iteration(db==k(s)).pka;
            norm = max(pka(1,:));
            plot([1 size(pka,2)], [0 0], ':', 'color', 'k')
            hold on;
            plot(1:size(pka,2), pka(1,:)/norm, '-', 'color', 'k')
            hold on;
            plot(1:size(pka,2), pka(3,:)/norm, '-', 'color', green)
            hold on;
            plot(1:size(pka,2), pka(2,:)/norm, '-', 'color', yellow)
            title(k(s)/noise)
            xlim([1 size(pka,2)])
            ylim([-0.5 2])
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        end
end
colormap(jet)
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
set(h, 'Name', model, 'NumberTitle', 'off')

function [delta, h, l] = normPKAd(pka, binoption)
pka0 = pka(1,:); % average
pka1 = pka(2,:); % high confidence
pka2 = pka(3,:); % low confidence
temp = (pka1 - pka2)/max(pka0); 
% binning
if binoption == 1
    temp = tcbin(temp, 4);
    pka0 = tcbin(pka0, 4);
    pka1 = tcbin(pka1, 4);
    pka2 = tcbin(pka2, 4);
end
delta = temp(end);
h = pka1(end)/max(pka0);
l = pka2(end)/max(pka0);