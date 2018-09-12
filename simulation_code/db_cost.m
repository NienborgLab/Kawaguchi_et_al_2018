function para = db_cost
%%
% estimate the cost of having decision bound on performance
hdx = 0.3*[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
kernel = [0.1 0.5 0.2 0.05 0 -0.05 -0.2 -0.5 -0.1];
len_tr = 5000;
len_frame = 100;
stmstrength = [0 3.125 6.25 12.5 25 50];
noise = 0.228;

% seed
rng(19891220);
    
% generate decision variables
hdxmat = nan(len_tr*2*length(stmstrength), len_frame);
sig = nan(len_tr*2*length(stmstrength),1); 
begin = 1;
for s = 1:length(stmstrength)
    % generate artificial dynamic stimuli
    restper = 100 - stmstrength(s);
    pd1 = (restper/length(hdx))*ones(length(hdx),1);
    pd2 = pd1;
    pd1(3) = stmstrength(s) + pd1(3);
    pd2(7) = stmstrength(s) + pd2(7);
    
    % estimate decision variables in a signal detection theory way
    stm = datasample(hdx, len_tr*len_frame, 'Replace', true, 'Weights', pd1);
    stm1 = reshape(stm, [len_tr, len_frame]); 
    stm = datasample(hdx, len_tr*len_frame, 'Replace', true, 'Weights', pd2);
    stm2 = reshape(stm, [len_tr, len_frame]);     
    hdxmat(begin:begin+len_tr*2-1,:) = [stm1; stm2];

    % signal    
    sig(begin:begin+len_tr*2-1,:) = [-stmstrength(s)*ones(len_tr,1); stmstrength(s)*ones(len_tr,1)];    
    
    begin = begin+len_tr*2;
end
% stm to dv
idv = hdxmat;
for s = 1:length(hdx)
    idv(stm==hdx(s)) = kernel(s);
end

% add noise
idv = hdxmat + normrnd(0, noise, size(hdxmat));

% integrated DV
dv = cumsum(idv, 2);
disp('DV generated')

% noise on stimulus
noisestm = std(hdxmat(:));
noiseidv = std(idv(:));

% psychometric functions & expected rewards
db = [[0:20]*noiseidv, inf];
% close all;
% h = figure;
col = colormap(jet(length(db)));
unistm = 100*unique(sig);
label = cell(1, length(db));
hitr = nan(1, length(db));
for d = 1:length(db)
    dv_temp = dv;
    ch = zeros(len_tr*2*length(stmstrength),1);
    if db(d) > 0
        dt = len_frame*ones(len_tr*2*length(stmstrength), 1);
        for l = 1:len_tr*2*length(stmstrength)
            dt_temp = find(abs(dv(l,:)) >= db(d), 1, 'first');
            if isempty(dt_temp)
                ch(l) = sign(dv_temp(l, end));
            else
                ch(l) = sign(dv_temp(l, dt_temp));
                dt(l) = dt_temp;
            end
        end
    else
        dt = zeros(len_tr*2*length(stmstrength), 1);
    end
    ch(ch==0) = datasample([-1 1], sum(ch==0));
    [pm, acc, mout] = compute_pm(sig, ch);
    label{d} = {['db:' num2str(db(d))]};    
    hitr(d) = 100*sum(dt < len_frame)/length(dt);
    disp('-------------------------------------------')
    disp(['uncorrected db: ' num2str(db(d)) ...
        ': bound hit: ' num2str(100*sum(dt < len_frame)/length(dt)) ' (%)' ...
        ', threshold: ' num2str(mout.params(2))])
    
    para.db(d).db = db(d);
    para.db(d).pm = pm;
    para.db(d).acc = acc;
    para.db(d).pm_fit = mout;
    
    if db(d) > 1000
        db(d) = 22*noiseidv;
    end
    db(d) = db(d)/noiseidv;
    
    % PM
    subplot(1,2,1)
    plot(100*mout.fitx/noisestm, mout.fity, '-', 'color', col(d,:), 'linewidth', 0.5)
    hold on;
    plot(unistm/noisestm, pm, 'o', 'color', col(d,:),'linewidth',0.5);    
    hold on;
    
    % expected rewards (percent correct)
    subplot(1,2,2)
    scatter(db(d), mean(acc), 50, 'o', 'markerfacecolor', col(d,:), ...
        'markeredgecolor','w');
    hold on;    
end
% shiny
subplot(1,2,1)
xx = [unistm(1) unistm(end)]/noisestm;
plot(xx, 0.5*ones(1,2), ':k')
hold on;
plot([0 0], [0 1], ':k')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
xlabel('% signal')
ylabel('probability of choice 1')
xlim(xx)
set(gca, 'XTick', [xx(1) 0 xx(2)], ...
    'XTickLabel', [-stmstrength(end) 0 stmstrength(end)])
set(gca, 'YTick', [0 0.5 1])

subplot(1,2,2)
yy = [0.48 0.8];
% decision bound range obtained from the fitting
fill([8.2385 8.769 8.769 8.2385 8.2385], [yy(1) yy(1) yy(2) yy(2) yy(1)], 'r', 'facealpha', 0.2, 'edgecolor', 'r', 'edgealpha', 0.1)
hold on;
ylim(yy)
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
xlabel('decision bound')
ylabel('probability correct')
set(gca, 'XTick', round([db(1) db(end-1) db(end)]), 'XTickLabel', round([db(1) db(end-1) inf]))
set(gca, 'YTick', [0.5 yy(2)])

set(gcf, 'Name', 'decision bound and performance', 'NumberTitle', 'off')

% subfunction
function [pm, acc, mout] = compute_pm(sig, ch)
unistm = unique(sig);
pm = 0.5*ones(1, length(unistm));
acc = pm;
n_resp = pm;
for i = 1:length(unistm)
    n_resp(i) = sum(sig==unistm(i));
    if abs(unistm(i)) > 0
        ch_temp = ch(sig==unistm(i));
        pm(i) = sum(sign(ch_temp)==1)/n_resp(i);
        acc(i) = sum(sign(ch_temp)==sign(unistm(i)))/n_resp(i);
    end    
end
% fit PM with a cumulative Gaussian
mout = fitPM(unistm', pm, n_resp,'Gaussian','MLE',0);