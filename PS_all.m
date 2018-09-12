function PS_all(focus)
% reproduce main results in Kawaguchi et al., 2018
%
% INPUT:
% focus ... 'pm','tc','learning', 'confsig','pka', 'tc_control', 'analysiswindow', or 'all'
%
% EXAMPLE: PS_all('all')
%

% path
mypath = cd;
addpath(genpath(mypath))

option = 'paper';
saveflag = 1; % default: save figures

% dealing with input
if nargin < 1; focus = 'all'; end

% yellow and green
yellow = [0.9576    0.7285    0.2285];
green = [0.1059    0.4706    0.2157];

% blue and red
red = [0.7922 0 0.1255];
blue = [0.0196 0.4431 0.6902];
purple = (red + blue)/2;

clc
close all;

% loop for animals
animals = {'mango','kiwi'};
animalsp = {'Animal B','Animal A'};
% signal strength in each animal
dcs = {[0    0.0300    0.0600    0.1250    0.2500    0.5000    1.0000],...
    [0    0.0310    0.0625    0.1250    0.2500    0.5000]};
for a = 1:length(animals)
    %%
    % load data
    disp([animals{a} '----------------'])
    disp('loading data......wait......')
    try
        load([mypath '/data/trmat_bandpass_' animals{a} '.mat'])       
    catch
        error('No data found. Please download the data by following the instruction in the https://github.com/NienborgLab/Kawaguchi_et_al_2018')
    end
    nses = length(unique(trmat(:,1)));
    disp(['The number of sessions: ' ...
        num2str(nses)])
    ntr = size(trmat,1);
    disp(['The number of trials: ' num2str(ntr)])

    %%
    % psychometric function (Figure 2b)
    if strcmp(focus, 'all') || strcmp(focus, 'pm')
       figure(1);
       subplot(1,2,a) 
       plot([-100 100], 0.5*[1 1], ':k')
       hold on;
       plot([0 0], [0 1], ':k')
       hold on;
       x = unique(trmat(:,3));
       lenx = length(x);
       ymean = zeros(1, lenx); 
       ysem = ymean;
       n = ysem;
       for s = 1:lenx
           n(s) = sum(trmat(trmat(:,3)==x(s),5)==1);
           ymean(s) = mean(trmat(trmat(:,3)==x(s),5), 1);
           ysem(s) = std(trmat(trmat(:,3)==x(s),5), [], 1)/sqrt(n(s));
       end
       n(x==-1 | x==1) = 0; 
       ymean(x==-1 | x==1) = nan;
       ysem(x==-1 | x==1) = nan;
       pmfit = fitPM(100*x',ymean,n,'Gaussian','MLE',0);             
       hold on;
       plot(pmfit.fitx, pmfit.fity, '-k')
       hold on;
       errorbar(100*x', ymean, ysem, 'k', 'linestyle', 'none', 'CapSize', 0)
       hold on;
       scatter(100*x', ymean, 15, 'o', 'markerfacecolor', ...
           'w', 'markeredgecolor', 'k')
       xlabel('signal strength (%)')
       ylabel('P(far choice)')
       ylim([0 1])
       title({animals{a}, [num2str(ntr) ' trials'],...
           ['threshold: ' num2str(pmfit.params(2))]})
       set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
       
       if a==2 && saveflag==1  
           savefig(gcf, [mypath '/data/figure_storage/PMs_' option '.fig'])
       end
    end     
    
    %%
    % signatures of confidence (Figure 5c, d)
    if strcmp(focus,'all') || strcmp(focus,'confsig')
        unises = unique(trmat(:,1));
        conflab = {'low', 'high'};
        acclab = {'error', 'correct'};
        if a == 1 && strcmp(option, 'paper')
            unises = unises(end-39:end);
        end
        nses = length(unises);
        binsize = 12;
        for n = 1:nses
            sesmat = trmat(trmat(:,1)==unises(n),:);
            for v = 1:2
                tempmat = sesmat(sesmat(:,9)==v-1,:);
                mps = tempmat(:,7);
                acc = tempmat(:,6);
                stm = tempmat(:,3);
                signature = confsig_compute(mps, acc, stm, dcs{a});
                confsig.animal(a).avrew(v).dc = dcs{a};
                confsig.animal(a).avrew(v).ntr(n,1) = signature.ntr;                    
                confsig.animal(a).avrew(v).cfx(n,:) = signature.cfx;
                confsig.animal(a).avrew(v).acx(n,:) = signature.acx;                
                for k = 1:2
                    confsig.animal(a).avrew(v).(['pm_' conflab{k} 'conf'])(n,:)...
                        = signature.pm_conf(k,:);
                    confsig.animal(a).avrew(v).([acclab{k} '_conf'])(n,:)...
                        = signature.correctness_conf(k,:);
                end
            end
        end

        for v = 1:2
            % stats
            confsig = stats_confsig(confsig, a, v);

            figure(4+v);
            % signature 1
            subplot(2,3,1+3*(a-1))
            [me, sem] = weighted(confsig.animal(a).avrew(v).acx, confsig.animal(a).avrew(v).ntr);
            fill_between(1:binsize, me-sem, me+sem, [0 0 0], 0.5)
            hold on;
            plot(1:binsize, me, '-k', 'linewidth',2)
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            ylabel('accuracy (%)')

            % signature 2
            subplot(2,3,2+3*(a-1))
            colors = {green, yellow};
            for k = 1:2
                [me, sem] = weighted(...
                    confsig.animal(a).avrew(v).(['pm_' conflab{k} 'conf']), confsig.animal(a).avrew(v).ntr);
                errorbar(1:length(dcs{a}), me, sem, '-','color',colors{k}, 'linewidth',0.75, 'CapSize', 0)
                hold on;
            end
            xlim([0.5 length(dcs{a})+0.5])
            set(gca, 'XTick', 1:length(dcs{a}), 'XTickLabel', dcs{a}) 
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

            % signature 3
            subplot(2,3,3+3*(a-1))
            colors = {0.4*[1 1 1], [0 0 0]};
            for k = 1:2
                [me, sem] = weighted(...
                    confsig.animal(a).avrew(v).([acclab{k} '_conf']), confsig.animal(a).avrew(v).ntr);
                errorbar(1:length(dcs{a}), me, sem, '-','color',colors{k}, 'linewidth',0.75, 'CapSize', 0);
                hold on;
            end
            xlim([0.5 length(dcs{a})+0.5])
            set(gca, 'XTick', 1:length(dcs{a}), 'XTickLabel', dcs{a})
            ylabel('confidence')
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            if a==2
                subplot(2,3,1+3*(a-1))
                xlabel('confidence')
                subplot(2,3,2+3*(a-1))
                xlabel('signal strength')
                subplot(2,3,3+3*(a-1))
                xlabel('signal strength')
            end
            if a==2
                set(gcf, 'Name', ['avrew. ' num2str(v-1)], 'NumberTitle', 'off')
                if saveflag == 1
                    savefig(gcf, [mypath '/data/figure_storage/confsig_avrew' num2str(v-1) '_' option '.fig'])
                end
           end
        end   
    end
    
    %% 
    % PKA (psychometric kernel amplitude) (Figure 2c, Figure 7a-d; datapoints)
    if strcmp(focus, 'all') || strcmp(focus, 'pka')    
        if a==1
            pka_all  = cell(1,3); pka_hc  = cell(1,3); pka_lc  = cell(1,3); 
        end
        pkantr = zeros(1, 2);
        unises = unique(trmat(:,1));
        % session groups due to animals' learning or experimental setups
        if a == 1
            sg{1} = 85-40:84;
        elseif a == 2
            sg{1} = unises(1:2);
            sg{2} = unises(3:38);
            sg{3} = unises(39:171);
            sg{4} = unises(172:213);
        end
        nbin = 4;
        pkmethod = 0;
        repeat = 1000; 
        nsg = length(sg);
        pka_all{a}  = zeros(2, nbin); pka_hc{a}  = zeros(2, nbin); pka_lc{a}  = zeros(2, nbin); 
        ntrsum = 0;
        for k = 1:nsg
            for v = 1:2
                trs = ismember(trmat(:,1), sg{k}) & trmat(:,9)==v-1 & trmat(:,3)==0;
                ntrs = sum(trs);
                ntrsum = ntrsum + ntrs;
                [pka_all_temp, pka_hc_temp, pka_lc_temp] = getPKA(trmat(trs,3), trmat(trs,end-450+1:end), ...
                    trmat(trs,5), trmat(trs,2), nbin, pkmethod, repeat);
                pka_all{a}(1,:) = pka_all{a}(1,:) + ntrs*pka_all_temp(1,:);
                pka_hc{a}(1,:) = pka_hc{a}(1,:) + ntrs*pka_hc_temp(1,:);
                pka_lc{a}(1,:) = pka_lc{a}(1,:) + ntrs*pka_lc_temp(1,:);
                if repeat > 0
                    pka_all{a}(2,:) = pka_all{a}(2,:) + ntrs*pka_all_temp(2,:);
                    pka_hc{a}(2,:) = pka_hc{a}(2,:) + ntrs*pka_hc_temp(2,:);
                    pka_lc{a}(2,:) = pka_lc{a}(2,:) + ntrs*pka_lc_temp(2,:);
                end
            end
        end
        pka_all{a} = pka_all{a}./ntrsum;
        pka_hc{a} = pka_hc{a}./ntrsum;
        pka_lc{a} = pka_lc{a}./ntrsum;
        pkantr(a) = ntrsum;
        
        % visualize
        figure(7);
        loc = [2 1];
        subplot(1,3, loc(a))
        hold on;
        norm = max(pka_all{a}(1,:));
        errorbar(1:nbin, pka_all{a}(1,:)/norm, pka_all{a}(2,:)/norm, '-k', 'CapSize', 0)
        hold on;
        errorbar(1:nbin, pka_lc{a}(1,:)/norm, pka_lc{a}(2,:)/norm, '-', 'color', green, 'CapSize', 0)
        hold on;
        errorbar(1:nbin, pka_hc{a}(1,:)/norm, pka_hc{a}(2,:)/norm, '-', 'color', yellow, 'CapSize', 0)
        xlim([0.5 nbin+0.5])
        title(animalsp{a})
        set(gca, 'XTick', [1 nbin], 'XTickLabel', [0 1.5]) 
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        if a==2
            ylabel('normalized PKA')
            for k = 1:2
                pka_all{end}(k,:) = (pkantr(1)*pka_all{1}(k,:) + pkantr(2)*pka_all{2}(k,:))/sum(pkantr);
                pka_lc{end}(k,:) = (pkantr(1)*pka_lc{1}(k,:) + pkantr(2)*pka_lc{2}(k,:))/sum(pkantr);
                pka_hc{end}(k,:) = (pkantr(1)*pka_hc{1}(k,:) + pkantr(2)*pka_hc{2}(k,:))/sum(pkantr);
            end
            subplot(1,3,3)
            hold on;
            norm = max(pka_all{end}(1,:));
            errorbar(1:nbin, pka_all{end}(1,:)/norm, pka_all{end}(2,:)/norm, '-k', 'CapSize', 0)
            hold on;
            errorbar(1:nbin, pka_lc{end}(1,:)/norm, pka_lc{end}(2,:)/norm, '-', 'color', green, 'CapSize', 0)
            hold on;
            errorbar(1:nbin, pka_hc{end}(1,:)/norm, pka_hc{end}(2,:)/norm, '-', 'color', yellow, 'CapSize', 0)
            xlim([0.5 nbin+0.5])
            title('both')
            set(gca, 'XTick', [1 nbin], 'XTickLabel', [0 1.5]) 
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        elseif a==1
            xlabel('time after stimulus onset (sec)')
        end

        if a==2 && saveflag==1
            set(gcf, 'Name', ['PKA_' option], 'NumberTitle', 'off')
            savefig([mypath '/data/figure_storage/PKA_' option '.fig'])
        end
    end
    
   %%
    % pupil size time-course (Figure 3a-c)
    if strcmp(focus,'all') || strcmp(focus,'tc')
        % load data
        disp([animals{a} '----------------'])
        disp('loading pupil data......wait......')
        
        % lowpass
        load([mypath '/data/psmat_lowpass_' animals{a} '.mat'])

        % quintile
        sesIDs = unique(psmat(:,1));
        quin = ones(size(psmat,1), 1);
        for n = 1:length(sesIDs)
            sesmat = psmat(psmat(:,1)==sesIDs(n) & psmat(:,3)==0 & psmat(:,9)==0, :);
            qui = prctile(sesmat(:,10), [0 20 40 60 80 100]);
            for t = 1:5
                quin(psmat(:,1)==sesIDs(n) & ...
                    psmat(:,10) >= qui(t) & qui(t+1) > psmat(:,10)) = t;
            end
        end
        
        % time-during the session
        figure(2000);
        loc = [2,1];
        subplot(4,2,loc(a))
        col = copper(5);
        colors =  cell(1,5);
        for t = 1:5
            colors{t} = col(t,:);
        end
        pstc_plot(psmat, {quin==1, quin==2, quin==3, quin==4, quin==5}, ...
            colors, {'-','-','-','-','-'})
        title(animalsp{a})

        % bandpass
        load([mypath '/data/psmat_bandpass_' animals{a} '.mat'])
        
        % available reward size
        subplot(4,2,loc(a)+2)
        cond0 = psmat(:,3) == 0 & psmat(:,9)==0 & psmat(:,8)==0;
        cond1 = psmat(:,3) == 0 & psmat(:,9)==1 & psmat(:,8)==0;
        pstc_plot(psmat, {cond0, cond1}, {blue, red},{'-','-'})
        ylabel('z-scored pupil size')

        % available reward size (nsplit = 3)
        subplot(4,2,loc(a)+4)
        cond0 = psmat(:,3) == 0 & psmat(:,4)==1 & psmat(:,8)==0;
        cond1 = psmat(:,3) == 0 & psmat(:,4)==0 & psmat(:,8)==0;
        cond2 = psmat(:,3) == 0 & psmat(:,4)==2 & psmat(:,8)==0;
        pstc_plot(psmat, {cond0, cond1, cond2}, {purple, blue, red},{'-','-','-'})
        ylabel('z-scored pupil size')

        % signal strength
        subplot(4,2,loc(a)+6)
        cond0 = psmat(:,11)==0 & psmat(:,9)==0 & psmat(:,8)==0;
        cond1 = psmat(:,11)==2 & psmat(:,9)==0 & psmat(:,8)==0;
        pstc_plot(psmat, {cond0, cond1}, {green, yellow},{'-','-'})
        xlabel('time after stimulus onset (ms)')

        set(gcf, 'Name', 'ps time-course', 'NumberTitle', 'off')

        if a==2 && saveflag==1
           savefig(gcf, [mypath '/data/figure_storage/PSTC_' option '.fig'])
        end

        clearvars psmat
    end
    
   %%
    % learning and pupil size (Figure 4a,b)
    if strcmp(focus, 'all') || strcmp(focus, 'learning')
        % load data
        disp([animals{a} '----------------'])
        disp('loading pupil data......wait......')
        if a < 3
            load([mypath '/data/psmat_bandpass_' animals{a} '.mat'])
        else
            continue
        end
        unises = unique(psmat(:,1));
        nses = length(unises);
        offset = 30;
        rocmat = nan(nses, 450+offset*2);
        rocmps = nan(nses, 1);
        threshold = nan(nses, 1);
        for n  = 1:nses
           % psychophysical threshold
           x = unique(psmat(psmat(:,1)==unises(n), 3));
           lenx = length(x);
           ymean = zeros(1, lenx); 
           k = ymean;
           for s = 1:lenx
               k(s) = sum(psmat(psmat(:,1)==unises(n) & psmat(:,3)==x(s), 5)==1);
               ymean(s) = mean(psmat(psmat(:,1)==unises(n) & psmat(:,3)==x(s), 5), 1);
           end
           k(x==-1 | x==1) = 0; 
           ymean(x==-1 | x==1) = nan;
           ysem(x==-1 | x==1) = nan;
           pmfit = fitPM(100*x',ymean,k,'Gaussian','MLE',0);        
           threshold(n) = pmfit.params(2);

           % area under ROC: easy vs hard in pupil size
           cond0 = psmat(:,1)==unises(n) & psmat(:,11)==0 & psmat(:,9)==0 & psmat(:,8)==0;
           cond1 = psmat(:,1)==unises(n) & psmat(:,11)==2 & psmat(:,9)==0 & psmat(:,8)==0;
           hardmat = psmat(cond0, end-540-offset+1:end-90+offset);
           easymat = psmat(cond1, end-540-offset+1:end-90+offset);
           parfor i = 1:size(hardmat, 2)
                rocmat(n,i) = rocN(easymat(:,i), hardmat(:,i));
           end     
           rocmps(n) = nanmean(rocmat(n, end - offset - 75 + 1 : end - offset)); 

            % visualization
            figure(1234);
            % ROC
            loc = [2,1];
            subplot(2,4,[1 2]+4*(loc(a)-1))
            imagesc(([1:450+2*offset]-offset)/300, 1:nses, rocmat)
            colormap(pink)
            hold on;
            c = colorbar('eastoutside');
            c.Label.String = 'aROC (easy vs hard trials)';
            yy = get(gca, 'YLim');
            plot([0 0], yy, '-k')
            hold on;
            plot(1.5*[1 1], yy, '-k')
            xlabel('time from stimulus onset (ms)')
            ylabel('session')
            title(animalsp{a})
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

            % threshold
            subplot(2,4,3+4*(loc(a)-1))
            plot(threshold, 1:nses, 'ok', 'markersize', 2)
            xlabel('psychophysical threshold (%)')
            ylabel('session')
            ylim([1 nses])
            set(gca, 'YDir', 'reverse')
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            
            % roc for mps vs threshold
            subplot(2,4,4+4*(loc(a)-1))
            scatter(rocmps, threshold, 10, 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'k')
            hold on;
            yy = get(gca, 'YLim');
            plot(0.5*[1 1], yy, '--k', 'linewidth', 0.5)
            beta = glmfit(rocmps, threshold, 'normal', 'link', 'identity', 'constant', 'on');
            x = min(rocmps):0.01:max(rocmps);
            hold on;
            plot(x, glmval(beta, x, 'identity'), '-k')
            xlabel('aROC in mean pupil size (easy vs hard)')
            ylabel('psychophysical threshold (%)')
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')            
                
            % autosave
            if a==2
                set(gcf, 'Name', 'learning', 'NumberTitle', 'off')
                if saveflag==1
                    savefig(gcf, [mypath '/data/figure_storage/PSvsLearning_' option '.fig'])
                end
            end
            clearvars psmat
        end
    end
        
    %%
    % control analysis on analysis window (Figure 6)
    if (strcmp(focus, 'all') || strcmp(focus, 'analysiswindow')) && strcmp(option, 'paper')
        % load data
        disp([animals{a} '----------------'])
        disp('loading pupil data......wait......')
        load([mypath '/data/psmat_bandpass_' animals{a} '.mat'])

        % analysis window
        samprate = 300;
        fixdur = 1.5;
        windowsize = [50, 150, 250, 350, 450]; % ms
        lenw = length(windowsize);
        stepf = 50; % ms
        offset = 0;
        conf0 = psmat(psmat(:,11) == 0 & psmat(:,9)==0 & psmat(:,8)==0, ...
            end-540-offset+1:end-90+offset);
        conf1 = psmat(psmat(:,11) == 2 & psmat(:,9)==0 & psmat(:,8)==0, ...
            end-540-offset+1:end-90+offset);  

        % iteration
        tc_cntr.animal(a).size = windowsize;
        tc_cntr.animal(a).ttest2_pval = nan(lenw, samprate*fixdur);
        tc_cntr.animal(a).stats1_spearman_p = nan(lenw, samprate*fixdur);
        tc_cntr.animal(a).stats2_permutation = nan(lenw, samprate*fixdur);
        tc_cntr.animal(a).stats3_spearman_p = nan(lenw, samprate*fixdur);
        for w = 1:lenw
            begin = 1;
            inc = (windowsize(w)*samprate*0.001):(stepf*samprate*0.001):450;
            for s = 1:length(inc)
                % p-val
                [~,p] = ttest2(nanmean(conf0(:, begin:begin+windowsize(w)*samprate*0.001 -1), 2), ...
                    nanmean(conf1(:, begin:begin+windowsize(w)*samprate*0.001 -1), 2));
                tc_cntr.animal(a).ttest2_pval(w, begin:begin+windowsize(w)*samprate*0.001-1) = p;
                
                % confidence signatures & its stats
                unises = unique(psmat(:,1));
                if a == 1 && strcmp(option, 'paper')
                    unises = unises(end-39:end);
                end
                nses = length(unises);
                for n = 1:nses
                    sesmat = psmat(psmat(:,1)==unises(n) & psmat(:,9)==0,:);
                    tcmat = sesmat(:, end-540-offset+1:end-90+offset);
                    mps = nanmean(tcmat(:, begin:begin+windowsize(w)*samprate*0.001 -1), 2);
                    acc = sesmat(:,6);
                    stm = sesmat(:,3);
                    signature = confsig_compute(mps, acc, stm, dcs{a});
                    confsig_aw.animal(a).avrew(1).dc = dcs{a};
                    confsig_aw.animal(a).avrew(1).ntr(n,1) = signature.ntr;                    
                    confsig_aw.animal(a).avrew(1).cfx(n,:) = signature.cfx;
                    confsig_aw.animal(a).avrew(1).acx(n,:) = signature.acx;
                    conflab = {'low', 'high'};
                    acclab = {'error', 'correct'};
                    for k = 1:2
                        confsig_aw.animal(a).avrew(1).(['pm_' conflab{k} 'conf'])(n,:)...
                            = signature.pm_conf(k,:);
                        confsig_aw.animal(a).avrew(1).([acclab{k} '_conf'])(n,:)...
                            = signature.correctness_conf(k,:);
                    end                    
                end
                confsig_aw = stats_confsig(confsig_aw, a, 1);
                tc_cntr.animal(a).stats1_spearman_p(w, begin:begin+windowsize(w)*samprate*0.001-1)...
                    = confsig_aw.animal(a).avrew(1).stats1_spearman(2);
                tc_cntr.animal(a).stats2_permutation(w, begin:begin+windowsize(w)*samprate*0.001-1)...
                    = confsig_aw.animal(a).avrew(1).stats2_permutation;
                tc_cntr.animal(a).stats3_spearman_p(w, begin:begin+windowsize(w)*samprate*0.001-1)...
                    = confsig_aw.animal(a).avrew(1).stats3_spearman(2);
                
                begin = begin + stepf*samprate*0.001;
            end
            tc_cntr.animal(a).window(w).confsig = confsig_aw.animal(a).avrew(1);
        end
        
        clearvars psmat
        
        % visualization
        if a==2
            figure(300);
            plot_analysiswindow_stats(tc_cntr);
            if saveflag==1
                savefig(gcf, [mypath '/data/figure_storage/analysiswindow.fig'])
            end
        end
    end   
    
    clearvars trmat
end

% subfunctions
function pstc_plot(psmat, conditions, colors, linetypes)
% data extraction
lenc = length(conditions);
nframe = 450;
sr = 300;
offset = 60;
unises = unique(psmat(:,1));
nses = length(unises);
tcmat = cell(1, lenc);
try
    for c = 1:lenc
        tcmat{c} = nan(nses, nframe+offset*2);
        for n = 1:nses
            tcmat{c}(n,:) = ...
                nanmean(psmat(conditions{c} & psmat(:,1)==unises(n), end-540-offset+1:end-90+offset), 1);
        end        
    end
    xx = [-0.2 1.7];
catch
    offset = 30;
    for c = 1:lenc
        tcmat{c} = nan(nses, nframe+offset*2);
        for n = 1:nses
            tcmat{c}(n,:) = ...
                nanmean(psmat(conditions{c} & psmat(:,1)==unises(n), end-540-offset+1:end-90+offset), 1);
        end        
    end
    xx = [-0.05 1.52];
end
% visualization
l = zeros(1,lenc);
for c = 1:lenc
    hold on;
    me = nanmean(tcmat{c}, 1);
    sem = nanstd(tcmat{c}, [], 1)/sqrt(nses);
    me = boxcar_smooth(me, 30);
    sem = boxcar_smooth(sem, 30);
    fill_between(([1:length(me)]-offset)./sr, me-sem, me+sem, colors{c}, 0.4)
    hold on;
    l(c) = plot(([1:length(me)]-offset)./sr,me, 'color',colors{c}, 'linestyle', linetypes{c}, 'linewidth', 0.5);
end
% stats
p_vals = nan(2, nframe+offset*2);
xlim(xx)
yy = get(gca, 'YLim');
if lenc > 1
    for t = 1:nframe
        statsmat = nan(nses, lenc);
        for c = 1:lenc
            statsmat(:,c) = tcmat{c}(:, t+offset);
        end
        if lenc > 2
            p_vals(1,t+offset) = anova1(statsmat, [], 'off');
            p_vals(2,t+offset) = kruskalwallis(statsmat, [], 'off');
        elseif lenc==2
            [~,p_vals(1,t+offset)] = ttest2(statsmat(:,1), statsmat(:,2));
            p_vals(2,t+offset) = ranksum(statsmat(:,1), statsmat(:,2));
        end
    end
    pvals = nan(size(p_vals));
    pvals(p_vals < 0.05/nframe) = 1;
    hold on;
    plot(([1:length(me)]-offset)./sr, (yy(1)+0.1*(yy(2)-yy(1)))*pvals(1,:), 'color', 0.8*ones(1,3), 'linewidth', 2);
end
hold on;
yy = get(gca, 'YLim');
yy(1) = floor(100*yy(1))/100;
yy(2) = ceil(100*yy(2))/100;
plot([0 0]./300, yy, '-k', 'linewidth', 0.5)
hold on;
plot([450 450]./300, yy, '-k', 'linewidth', 0.5)
ylim(yy)
set(gca, 'XTick', [0 1.5])
set(gca,'YTick', yy)
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

function newvec = boxcar_smooth(oldvec, carlen)
% apply boxcar convolution to smooth the given vector
% this function automatically removes edge noise

car = ones(1, carlen)./carlen;
n = ceil(carlen/2);
newvec = conv([oldvec(1:n) oldvec oldvec(end-n+1:end)], car, 'same');
newvec(1:n) = []; newvec(end-n+1:end) = [];

function signature = confsig_compute(mps, acc, stm, dcrange)
% adjust stimulus range to dc range
sign_dc = unique([-dcrange dcrange]);
dc = unique(stm);
lendc = length(dc);
for d = 1:lendc
    [~,idx] = min(abs(sign_dc - dc(d)));
    stm(stm==dc(d)) = sign_dc(idx);
end
dc = dcrange;
lendc = length(dcrange);

% signature 1
binsize = 12;
ntr = length(mps);
frameperbin = floor(ntr/binsize);
[sorted_cf, si] = sort(mps);
sorted_acc = acc(si);
cfx = zeros(1, binsize);
acx = cfx;
begin = 1;
for k = 1:binsize
    cfx(k) = nanmean(sorted_cf(begin:begin+frameperbin-1));
    vec = sorted_acc(begin:begin+frameperbin-1);
    acx(k) = nansum(vec==1)/length(vec);
    begin = begin + frameperbin;
end
signature.cfx = cfx;
signature.acx = acx;

% signature 2 & 3
loc = [-1, 1];
for s = 1:lendc
    pcor = zeros(2,2);
    ndc = zeros(1,2);
    confcor = zeros(2,2);
    for l = 1:2
        % signature 2
        med = nanmedian(mps(stm==loc(l)*dc(s)));
        ndc(l) = nansum(stm==loc(l)*dc(s));
        pcor(l,1) = nansum(mps < med & stm==loc(l)*dc(s) & acc==1) ...
            /nansum(mps < med & stm==loc(l)*dc(s));
        pcor(l,2) = nansum(mps > med & stm==loc(l)*dc(s) & acc==1) ...
            /nansum(mps > med & stm==loc(l)*dc(s));
        
        % signature 3
        confcor(l,1) = nanmean(mps(stm==loc(l)*dc(s) & acc==0));
        confcor(l,2) = nanmean(mps(stm==loc(l)*dc(s) & acc==1));
    end
    % signature 2 & 3
    if sum(isnan([pcor(1,:), confcor(1,:)])) > 0
        signature.pm_conf(1,s) = nanmean(pcor(:,1), 1);
        signature.pm_conf(2,s) = nanmean(pcor(:,2), 1);
        signature.correctness_conf(1,s) = nanmean(confcor(:,1), 1);
        signature.correctness_conf(2,s) = nanmean(confcor(:,2), 1);
    else
        signature.pm_conf(1,s) = (pcor(1,1) + pcor(2,2))/2;
        signature.pm_conf(2,s) = (pcor(1,2) + pcor(2,2))/2;
        signature.correctness_conf(1,s) = (confcor(1,1) + confcor(2,2))/2;
        signature.correctness_conf(2,s) = (confcor(1,2) + confcor(2,2))/2;
    end
end

signature.signal = dc;
signature.ntr = ntr;

function confsig = stats_confsig(confsig, a, v)
% signature 1
cfx = mean(confsig.animal(a).avrew(v).cfx, 1)';
acx = mean(confsig.animal(a).avrew(v).acx, 1)';
[rr, pp] = corrcoef(cfx, acx);
confsig.animal(a).avrew(v).stats1_pearson = [rr(1,2), pp(1,2)];
[rho,pval] = corr(cfx, acx, 'type', 'Spearman');
confsig.animal(a).avrew(v).stats1_spearman = [rho, pval];

% signature 2
repeat = 1000;
nses = size(confsig.animal(a).avrew(v).ntr, 1);
thrediff = zeros(1, nses);
for r = 1:repeat
    rses = datasample(1:nses, nses, 'Replace', true);
    fitted0 = fitPM(confsig.animal(a).avrew(v).dc, ...
        mean(confsig.animal(a).avrew(v).pm_lowconf(rses,:), 1), ...
        ones(1, length(confsig.animal(a).avrew(v).dc)), 'Gaussian', 'MLE', 0);
    fitted1 = fitPM(confsig.animal(a).avrew(v).dc, ...
        mean(confsig.animal(a).avrew(v).pm_highconf(rses,:), 1), ...
        ones(1, length(confsig.animal(a).avrew(v).dc)), 'Gaussian', 'MLE', 0);
    thrediff(r) = fitted0.params(2) - fitted1.params(2);
end
confsig.animal(a).avrew(v).stats2_permutation...
    = sum(thrediff <= 0)/repeat;

% signature 3 (compare with a prediction by a SDT) 
if a == 1
    ncor = 7;
    nerr = 4;
    if v==1
        pred = [0.7476 0.7554 0.7606 0.7751 0.7994 0.8517 0.9480,...
            0.7521 0.7378 0.7372 0.7287 0.7072 0.6794 0.6476];
    elseif v==2
        pred = [0.7465 0.7538 0.7666 0.7729 0.8004 0.8602 0.9533,...
            0.7492 0.7459 0.7355 0.7256 0.7096 0.6774 0.6131];
    end
elseif a == 2
    ncor = 6;
    nerr = 3;
    if v==1
        pred = [0.7485 0.7656 0.7730 0.7966 0.8471 0.9409,...
            0.7477 0.7352 0.7300 0.7097 0.6869 0.6437];
    elseif v==2
        pred = [0.7498 0.7628 0.7778 0.8085 0.8701 0.9692,...
            0.7501 0.7374 0.7287 0.7036 0.6635 0.6395];
    end
end
pred = pred(1:ncor+nerr);
data = [nanmean(confsig.animal(a).avrew(v).correct_conf(:, 1:ncor), 1), ...
    nanmean(confsig.animal(a).avrew(v).error_conf(:, 1:nerr), 1)];
[rr, pp] = corrcoef(pred', data');
confsig.animal(a).avrew(v).stats3_pearson = [rr(1,2), pp(1,2)];
[rho,pval] = corr(pred', data', 'type', 'Spearman');
confsig.animal(a).avrew(v).stats3_spearman = [rho, pval];

function cp = rocN(x,y)
N = 100;
[m n] = size(x);
x = reshape(x,1,m*n);
[m n] = size(y);
y = reshape(y,1,m*n);

zlo = min([min(x(:)) min(y(:))]);
zhi = max([max(x(:)) max(y(:))]);
z = linspace(zlo,zhi,N);
fa = zeros(1,N);	% allocate the vector
hit = zeros(1,N);
for i = 1:N
  fa(N-i+1) = sum(y > z(i));
  hit(N-i+1) = sum(x > z(i));
end
[m,ny] = size(y);
fa = fa/ny;
[m,nx] = size(x);
hit = hit/nx;
fa(1) = 0;
hit(1) = 0;
fa(N) = 1;
hit(N) = 1;
cp = trapz(fa,hit);
