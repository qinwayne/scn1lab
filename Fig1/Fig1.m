clc
clear
close all

% Define groups and colors
WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
all_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63','01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};

num_WT = length(WT_fish);
num_all = length(all_fish);

cb_blue = [0, 0, 255] / 255;
cb_red  = [255, 0, 0] / 255;

% Load data and parameters
load('seizure_time.mat')
win = 6;
ignsec = 2;

% Extract valid seizure times
seizuretime = cell(1, num_all);
for fi = 1:num_all
    skewnessVals = seizureSkewness{fi};
    ci = seizureTime_ci{fi};
    ci = ci(abs(skewnessVals) > 0.8);

    cimat = ci + (-win:0.5:win)';
    mcmat = seizureTime_mc{fi} + (-win:0.5:win)';
    temp = unique(round(intersect(cimat(:), mcmat(:))));

    tm = seizureTime_tm{fi};
    if ~isempty(tm)
        tmmat = round(tm) + (-win:win)';
        temp = intersect(temp(:), tmmat(:));
    end

    if ~isempty(temp)
        valid_idx = [true, diff(temp(:)') > ignsec];
        seizuretime{fi} = temp(valid_idx) + win;
    end
end

% Create the main figure layout
figure('Color', 'w', 'Position', [100, 100, 600, 900]);

% Panel C: Seizure Number Comparison (Top)
subplot(3, 1, 1);
hold on;

wtno = [];
mutno = [];
for fi = 1:num_all
    if isempty(seizuretime{fi}), continue; end
    szinx = sum(seizuretime{fi} < 2500 & seizuretime{fi} > 350);
    if szinx == 0, continue; end

    if fi <= num_WT
        wtno(end+1) = szinx;
    else
        mutno(end+1) = szinx;
    end
end

data1 = [wtno'; mutno'];
group1 = [repmat({'WT'}, length(wtno), 1); repmat({'MUT'}, length(mutno), 1)];

violinplot(data1, group1, 'ShowMean', true, 'ShowBox', false, 'GroupOrder', {'WT', 'MUT'}, 'ViolinColor', [cb_blue; cb_red], 'MarkerSize', 20, 'LineWidth', 2);

[h1, p1] = ttest2(wtno, mutno);
if h1
    yl = ylim;
    plot([1, 2], [yl(2)*0.95, yl(2)*0.95], 'k', 'LineWidth', 1.5);
    text(1.5, yl(2)*0.98, get_sig_star(p1), 'FontSize', 14, 'HorizontalAlignment', 'center');
end

ylabel('Seizures per animal', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1, 'Box', 'off');
ylim([0, 20]);
xlim([0.5, 2.5]);

% Panel 2: Inter Seizure Interval Comparison (Middle)
subplot(3, 1, 2);
hold on;

wtisi = [];
mutisi = [];
for fi = 1:num_all
    if length(seizuretime{fi}) < 2, continue; end
    
    isi = diff(seizuretime{fi});
    isi = isi(:)'; % FIX: Force isi to always be a row vector
    
    if fi <= num_WT
        wtisi = [wtisi, isi];
    else
        mutisi = [mutisi, isi];
    end
end

wt_data = wtisi(wtisi > ignsec);
wt_data = wt_data(~isoutlier(wt_data));

mut_data = mutisi(mutisi > ignsec);
mut_data = mut_data(~isoutlier(mut_data));

data2 = [wt_data'; mut_data'];
group2 = [repmat({'WT'}, length(wt_data), 1); repmat({'MUT'}, length(mut_data), 1)];

violinplot(data2, group2, 'ShowMean', true, 'ShowBox', false, 'GroupOrder', {'WT', 'MUT'}, 'ViolinColor', [cb_blue; cb_red], 'MarkerSize', 20, 'LineWidth', 2);

[h2, p2] = ttest2(wt_data, mut_data);
if h2
    yl = ylim;
    plot([1, 2], [yl(2)*0.95, yl(2)*0.95], 'k', 'LineWidth', 1.5);
    text(1.5, yl(2)*0.98, get_sig_star(p2), 'FontSize', 14, 'HorizontalAlignment', 'center');
end

ylabel('Inter seizure interval (s)', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1, 'Box', 'off');
xlim([0.5, 2.5]);

% Panel E: Seizure Number Per Stage Comparison (Bottom)
subplot(3, 1, 3);
hold on;

stages_ptz = 400:300:4200;
binstages = [stages_ptz / 2, 2150];
szhist = zeros(num_all, length(binstages) - 1);

for fi = 1:num_all
    if isempty(seizuretime{fi}), continue; end
    [N, ~] = histcounts(seizuretime{fi}, binstages);
    szhist(fi, :) = N;
end

wt_hist = szhist(1:num_WT, :);
mut_hist = szhist(num_WT+1:end, :);

wtszs = mean(wt_hist, 1);
wtsigs = var(wt_hist, 0, 1) / num_WT;

mutszs = mean(mut_hist, 1);
mutsigs = var(mut_hist, 0, 1) / (num_all - num_WT);

data3 = [wtszs; mutszs]';
errors3 = [sqrt(wtsigs); sqrt(mutsigs)]';

b = bar(data3, 'BarWidth', 1);
b(1).FaceColor = cb_blue;
b(2).FaceColor = cb_red;

numGroups = size(data3, 1);
numBars = size(data3, 2);
groupWidth = min(0.6, numBars / (numBars + 1.5));

for i = 1:numBars
    x = (1:numGroups) - groupWidth/2 + (2*i - 1) * groupWidth / (2*numBars);
    errorbar(x, data3(:, i), errors3(:, i), 'k.', 'LineWidth', 1);
end

labelList = cell(1, length(binstages) - 1);
for i = 1:length(binstages) - 1
    val1 = binstages(i) - 300;
    labelList{i} = mat2str(val1);

    all_stage_data = [wt_hist(:, i); mut_hist(:, i)];
    group_labels = [repmat({'WT'}, num_WT, 1); repmat({'Mutant'}, num_all - num_WT, 1)];

    [~, ~, stats] = anova1(all_stage_data, group_labels, 'off');
    posthoc_results = multcompare(stats, 'CType', 'dunn-sidak', 'Display', 'off');

    p_val = posthoc_results(6);
    if p_val < 0.05
        text(i - 0.1, max(data3(i,:)) + max(errors3(i,:)) + 0.1, get_sig_star(p_val), 'FontSize', 15);
    end
end

ylabel('Seizures per animal', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
xticks(1:2:length(labelList) + 1);
xticklabels(labelList(1:2:end));
set(gca, 'FontSize', 12, 'LineWidth', 1, 'Box', 'off');

%% Helper Function for Significance Stars
function star = get_sig_star(p)
    if p < 0.001
        star = '***';
    elseif p < 0.01
        star = '**';
    elseif p < 0.05
        star = '*';
    else
        star = '';
    end
end