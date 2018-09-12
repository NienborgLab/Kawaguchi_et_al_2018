function plot_analysiswindow_stats(tc_cntr)
% plot results from control analysis on analysis window

% extract data
fn = {'ttest2_pval','stats1_spearman_p','stats2_permutation','stats3_spearman_p'};
clab = {'easy vs hard', 'pupil size vs accuracy', ...
    'psychophysical threshold split by pupil size', 'pupil size vs choice correctness'};
lenfn = length(fn);
animals = {'Animal B', 'Animal A'};
loc = [2,1];
% close all;
% h = figure;
cmap = flipud(parula(100));
cmap(end,:) = 0.5*ones(1,3);
for a = 1:2
    for f = 1:lenfn
        subplot(lenfn, 2, loc(a)+2*(f-1))
        imagesc(1:450, tc_cntr.animal(a).size, log(tc_cntr.animal(a).(fn{f})));
        colormap(cmap)
        bonferonni = 1;
        lim = caxis;
        if f==1            
            bonferonni = 450;
            if a==2
                lim(1) = -700;
            elseif a==1
                lim(1) = -60;
            end
        end
        lim(2) = log(0.05/bonferonni);
        caxis(lim)
        hold on;
        plot([435 450 450 435 435], [200 200 300 300 200], '-','color', [0 1 0.4],'linewidth',1.5)
        c = colorbar('northoutside');
        set(c, 'Ticks', lim, 'TickLabels', [bonferonni/exp(abs(lim(1))) 0.05])
        if a==2
            c.Label.String = clab{f};
        end
        if f==1
            title(animals{a})
        end
        if a==2 && f==2
            ylabel('analysis window length (ms)')
        end
        if a==2 && f==lenfn
            xlabel('time from stimulus onset (ms)')
        end
        set(gca,'YDir','normal')
        set(gca, 'XTick', [1 450], 'XTickLabel', [0 1500])
        set(gca, 'YTick', tc_cntr.animal(a).size)
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    end
end