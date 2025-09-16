% This script processes and visualizes seizure data from wild-type (WT) and
% SCN1a mutant (SCN1) zebrafish. It identifies seizure events, analyzes
% seizure numbers, inter-seizure intervals, and seizure timing over a PTZ
% (Pentylenetetrazol) stimulation period.
%
% Requires:
% - A 'seizure_time.mat' file containing seizure data.
% - The 'violinplot' toolbox from the MATLAB File Exchange or GitHub.
%
% Author: [Your Name]
% Date: [Date]
% GitHub: [Your GitHub URL]

%% 1. INITIALIZATION & DATA SETUP
% Clear workspace, close figures, and define experimental groups.
clc;
clear all;
close all;

% Define fish IDs for each group.
WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};
all_fish = [WT_fish, SCN1_fish];

% Define parameters for seizure detection and analysis.
slide_window_size = 300; % Sliding window size in seconds for binning.
ignored_inter_seizure_time = 60; % Time in seconds to ignore between detected seizures.

% Create time bins for seizure counting.
stages_ptz(1) = 400; % Start time.
i = 2;
while 400 + slide_window_size * (i-1) <= 4200
    stages_ptz(i) = 400 + slide_window_size * (i-1);
    i = i + 1;
end
stages_ptz = stages_ptz / 2; % Adjust to match data timescale.
bin_stages = [stages_ptz, 2150]; % Finalize bin edges.

% Load seizure data from the .mat file.
load('seizure_time.mat');

%% 2. SEIZURE TIME CONSOLIDATION
% Loop through different parameter combinations for robust detection.
% This nested loop structure allows for testing different window sizes and
% inter-seizure intervals.
for time_window_half_width = 6 % Test half-width of time window for intersection.
    for ignored_inter_seizure_time = 2 % Test ignored seconds between seizures.

        % Process each fish's data to find unified seizure times.
        seizure_time = cell(1, length(all_fish));
        novidfish = false(length(all_fish), 1);

        for fi = 1:length(all_fish)
            % Extract data for the current fish.
            ci = seizureTime_ci{fi};
            skewness_vals = seizureSkewness{fi};
            mc = seizureTime_mc{fi};
            tm = seizureTime_tm{fi};

            % Filter seizure times based on skewness threshold.
            ci = ci(abs(skewness_vals) > 0.8);

            % Create time intervals for intersection.
            cimat = ci + [-time_window_half_width:0.5:time_window_half_width]';
            mcmat = mc + [-time_window_half_width:0.5:time_window_half_width]';

            % Find times common to both CI and MC methods.
            temp_seizure_times = unique(round(intersect(reshape(cimat, 1, []), reshape(mcmat, 1, []))));

            % Intersect with TM data if available.
            if ~isempty(tm)
                tm_mat = round(tm) + [-time_window_half_width:time_window_half_width]';
                temp_seizure_times = intersect(reshape(temp_seizure_times, 1, []), reshape(tm_mat, 1, []));
            else
                % Mark fish without video data.
                novidfish(fi) = true;
            end

            % Store consolidated seizure times, filtering out events that are too close.
            if ~isempty(temp_seizure_times)
                seizure_time{fi} = temp_seizure_times([true, diff(temp_seizure_times) > ignored_inter_seizure_time]) + time_window_half_width;
            end
        end

        %% 3. VISUALIZATION & STATISTICAL ANALYSIS
        % This section generates several figures for analysis.

        % Define a colorblind-safe palette for plots.
        cb_blue = [0, 0, 255] / 255;
        cb_red = [255, 0, 0] / 255;

        % --- Figure 1: Seizure Number Comparison (Violin Plot) ---
        figure('Color', 'w');
        hold on;

        % Count seizures per animal for WT and Mutant groups.
        wt_seizure_count = [];
        mut_seizure_count = [];
        for fi = 1:length(all_fish)
            % Filter seizures within a specific time range (350s to 2500s).
            sz_count = length(find(seizure_time{fi} < 2500 & seizure_time{fi} > 350));
            if fi <= length(WT_fish)
                wt_seizure_count = [wt_seizure_count, sz_count];
            else
                mut_seizure_count = [mut_seizure_count, sz_count];
            end
        end

        % Remove entries with zero seizures.
        wt_seizure_count(wt_seizure_count == 0) = [];
        mut_seizure_count(mut_seizure_count == 0) = [];

        % Prepare data for violin plot.
        group = [repmat({'WT'}, length(wt_seizure_count), 1); repmat({'MUT'}, length(mut_seizure_count), 1)];
        data = [wt_seizure_count'; mut_seizure_count'];

        % Generate violin plot.
        vp = violinplot(data, group, 'ShowMean', true, 'ShowBox', false, 'GroupOrder', {'WT', 'MUT'}, 'ViolinColor', [cb_blue; cb_red], 'MarkerSize', 20, 'LineWidth', 2);

        % Perform t-test for statistical significance.
        [h, p] = ttest2(wt_seizure_count, mut_seizure_count);
        sig_star = '';
        if p < 0.001
            sig_star = '***';
        elseif p < 0.01
            sig_star = '**';
        elseif p < 0.05
            sig_star = '*';
        end

        % Annotate the plot with significance star.
        if h
            yl = ylim;
            plot([1, 2], [yl(2) * 0.95, yl(2) * 0.95], 'k', 'LineWidth', 1.5);
            text(1.5, yl(2) * 0.98, sig_star, 'FontSize', 14, 'HorizontalAlignment', 'center');
        end

        % Format plot axes and appearance.
        ylabel('Seizures per animal', 'FontSize', 12);
        set(gca, 'FontSize', 12, 'LineWidth', 1, 'Box', 'off');
        ylim([0, 20]);
        xlim([0.5, 2.5]);
        set(gcf, 'Position', [848, 882, 560, 232]);

        % --- Figure 2: Seizure Number Per Stage (Bar Plot) ---
        figure;
        
        % Create histograms of seizure times for each fish.
        sz_hist = zeros(length(all_fish), length(bin_stages) - 1);
        for fi = 1:length(all_fish)
            h = histogram(seizure_time{fi}, bin_stages, 'Visible', 'off');
            sz_hist(fi, :) = h.Values;
        end

        % Calculate mean and standard error for each group per stage.
        wt_means = mean(sz_hist(1:length(WT_fish), :), 1);
        wt_stderrs = sqrt(var(sz_hist(1:length(WT_fish), :), 1) / length(WT_fish));

        mut_means = mean(sz_hist(length(WT_fish)+1:end, :), 1);
        mut_stderrs = sqrt(var(sz_hist(length(WT_fish)+1:end, :), 1) / length(SCN1_fish));

        % Combine data and errors for bar plot.
        data_to_plot = [wt_means; mut_means]';
        errors_to_plot = [wt_stderrs; mut_stderrs]';

        % Create bar plot.
        b = bar(data_to_plot, 'BarWidth', 1);
        hold on;

        % Add error bars to the plot.
        num_groups = size(data_to_plot, 1);
        num_bars = size(data_to_plot, 2);
        group_width = min(0.6, num_bars / (num_bars + 1.5));
        for i = 1:num_bars
            x = (1:num_groups) - group_width / 2 + (2 * i - 1) * group_width / (2 * num_bars);
            errorbar(x, data_to_plot(:, i), errors_to_plot(:, i), 'k.', 'LineWidth', 1);
        end
        
        % Perform statistical tests (ANOVA and post-hoc) for each time bin.
        for i = 1:size(bin_stages, 2) - 1
            wt_data = sz_hist(1:length(WT_fish), i);
            mut_data = sz_hist(length(WT_fish)+1:end, i);
            
            % Skip if one group is empty.
            if isempty(wt_data) || isempty(mut_data)
                continue;
            end
            
            all_data = [wt_data; mut_data];
            group_labels = [repmat({'WT'}, length(WT_fish), 1); repmat({'Mutant'}, length(SCN1_fish), 1)];
            
            % Perform one-way ANOVA.
            [p, ~, stats] = anova1(all_data, group_labels, 'off');
            
            % Perform Šidák-corrected post hoc test.
            posthoc_results = multcompare(stats, 'CType', 'dunn-sidak', 'Display', 'off');
            p_posthoc = posthoc_results(1, 6); % p-value for the WT vs Mutant comparison.
            
            % Add significance stars to the plot.
            sig_star = '';
            if p_posthoc < 0.05 && p_posthoc >= 0.01
                sig_star = '*';
            elseif p_posthoc < 0.01 && p_posthoc >= 0.001
                sig_star = '**';
            elseif p_posthoc < 0.001
                sig_star = '***';
            end
            if p_posthoc < 0.05
                text(i - 0.1, max(data_to_plot(i, :)) * 1.1, sig_star, "FontSize", 15);
            end
        end
        
        % Format plot axes and labels.
        ylabel('Seizures per animal');
        xticks(1:size(bin_stages, 2) - 1);
        xlabel('Time (s)');
        xticklabels(bin_stages(1:end-1));
        set(gcf, 'position', [848, 882, 560, 232]);

        % --- Figure 3: Correlation Distribution Skewness Over Time ---
        figure;
        hold on;
        
        % Aggregate skewness data for each group.
        mut_x = []; mut_y = [];
        wt_x = []; wt_y = [];
        for fi = 1:length(seizureTime_ci)
            ci = seizureTime_ci{fi};
            skewness_vals = seizureSkewness{fi};
            
            % Filter data points with skewness > 0.5.
            x = ci(abs(skewness_vals) > 0.5);
            y = skewness_vals(abs(skewness_vals) > 0.5);
            
            if fi > length(WT_fish)
                mut_x = [mut_x, x];
                mut_y = [mut_y, y];
            else
                wt_x = [wt_x, x];
                wt_y = [wt_y, y];
            end
        end
        
        % Plot and fit linear model for Mutant group.
        c_mut = fitlm(mut_x, mut_y);
        h_mut = plot(c_mut);
        scatter(mut_x, mut_y, 'rx');
        text(600, -2.5, num2str(c_mut.Coefficients.Estimate(2), '%.4f'), 'color', 'r');
        
        % Plot and fit linear model for WT group.
        c_wt = fitlm(wt_x, wt_y);
        h_wt = plot(c_wt);
        fitHandle_wt = findobj(h_wt, 'DisplayName', 'Fit');
        fitHandle_wt.Color = cb_blue;
        cbHandles_wt = findobj(h_wt, 'DisplayName', 'Confidence bounds');
        set(cbHandles_wt, 'Color', cb_blue);
        scatter(wt_x, wt_y, 'bx');
        text(700, 0, num2str(c_wt.Coefficients.Estimate(2), '%.4f'), 'color', 'b');
        
        % Format plot axes and labels.
        xlim([300, 2100]);
        ylim([-3, 0]);
        xlabel('Time (s)');
        ylabel('Correlation distribution skewness');
        legend('off');
        set(gcf, 'position', [848, 882, 560, 232]);
        title([]);

        % --- Figure 4: Inter-Seizure Interval (ISI) Comparison ---
        figure('Color', 'w');
        hold on;

        % Calculate ISI for each group.
        wt_isi = [];
        mut_isi = [];
        for fi = 1:length(all_fish)
            isi = diff(seizure_time{fi});
            if fi <= length(WT_fish)
                wt_isi = [wt_isi, isi];
            else
                mut_isi = [mut_isi, isi];
            end
        end

        % Filter out ISIs shorter than the ignored time and outliers.
        wt_isi_filtered = wt_isi(wt_isi > ignored_inter_seizure_time);
        mut_isi_filtered = mut_isi(mut_isi > ignored_inter_seizure_time);
        wt_isi_filtered = wt_isi_filtered(~isoutlier(wt_isi_filtered));
        mut_isi_filtered = mut_isi_filtered(~isoutlier(mut_isi_filtered));

        % Prepare data for violin plot.
        group = [repmat({'WT'}, length(wt_isi_filtered), 1); repmat({'MUT'}, length(mut_isi_filtered), 1)];
        data = [wt_isi_filtered'; mut_isi_filtered'];

        % Generate violin plot.
        vp = violinplot(data, group, 'ShowMean', true, 'ShowBox', false, 'GroupOrder', {'WT', 'MUT'}, 'ViolinColor', [cb_blue; cb_red], 'MarkerSize', 20, 'LineWidth', 2);

        % Perform t-test for statistical significance.
        [h, p] = ttest2(wt_isi_filtered, mut_isi_filtered);
        sig_star = '';
        if p < 0.001
            sig_star = '***';
        elseif p < 0.01
            sig_star = '**';
        elseif p < 0.05
            sig_star = '*';
        end

        % Annotate the plot with significance star.
        if h
            yl = ylim;
            plot([1, 2], [yl(2) * 0.95, yl(2) * 0.95], 'k', 'LineWidth', 1.5);
            text(1.5, yl(2) * 0.98, sig_star, 'FontSize', 14, 'HorizontalAlignment', 'center');
        end

        % Format plot axes and appearance.
        ylabel('Inter-seizure interval (s)', 'FontSize', 12);
        set(gca, 'FontSize', 12, 'LineWidth', 1, 'Box', 'off');
        xlim([0.5, 2.5]);
        set(gcf, 'position', [848, 882, 560, 232]);

        % --- Figure 5: First Seizure Latency Comparison ---
        figure('Color', 'w');
        hold on;
        
        wt_first_seizure = [];
        mut_first_seizure = [];
        for fi = 1:length(all_fish)
            % Find the first seizure time after a specific buffer (350s).
            seizure_times = seizure_time{fi};
            if fi <= length(WT_fish)
                % Special case for fish '14' (index 10) to handle a known outlier.
                if fi == 10
                    seizure_times(1) = [];
                end
                wt_first_seizure = [wt_first_seizure, min(seizure_times(seizure_times > 350))];
            else
                mut_first_seizure = [mut_first_seizure, min(seizure_times(seizure_times > 350))];
            end
        end
        
        % Plot individual data points and boxcharts.
        scatter(ones(size(wt_first_seizure)), wt_first_seizure, 'b.');
        boxchart(ones(size(wt_first_seizure)), wt_first_seizure, 'BoxFaceColor', 'b');
        scatter(ones(size(mut_first_seizure)) * 2, mut_first_seizure, 'r.');
        boxchart(ones(size(mut_first_seizure)) * 2, mut_first_seizure, 'BoxFaceColor', 'r', 'MarkerColor', 'r', 'MarkerStyle', '.');
        
        % Perform t-test and annotate significance.
        [h, p] = ttest2(wt_first_seizure, mut_first_seizure);
        sig_star = '';
        if p < 0.05 && p >= 0.01
            sig_star = '*';
        elseif p < 0.01 && p >= 0.001
            sig_star = '**';
        elseif p < 0.001
            sig_star = '***';
        end
        
        if h
            yl = ylim;
            plot([1, 2], [yl(2) * 0.95, yl(2) * 0.95], 'k');
            text(1.5, yl(2) * 0.98, sig_star, 'FontSize', 14, 'HorizontalAlignment', 'center');
        end
        
        % Format plot axes and labels.
        xticks([1, 2]);
        xticklabels({'WT', 'SCN1'});
        ylabel('First seizure time (s)');
        set(gcf, 'position', [848, 882, 560, 232]);

    end
end