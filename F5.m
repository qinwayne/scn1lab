%% Panel 1
%
% This script statistically analyzes and visualizes differences in graph
% theory metrics between WT and mutant fish over multiple experimental stages.
% It performs ANOVA and then creates a bubble plot to show the
% significance and magnitude of the differences at each stage.


%% 1. INITIALIZATION AND DATA PREPARATION
clear all;
clc;
close all;

% Load the data file.
load('whole_GT24.mat');

% Ensure the field order is consistent.
orderedFields = fieldnames(gt);
gt = orderfields(gt, orderedFields);

% Manually create the 'gtom' matrix by averaging data from 'gtommap'.
% The `gtommap` cell array is assumed to contain a single metric value
% at the fourth column of each cell.
gt.gtom = zeros(35, 24);
for i = 1:35
    for j = 1:24
        gt.gtom(i, j) = nanmean(gt.gtommap{i, j}(:, 4));
    end
end

% Define the full list of fish and their genotypes.
WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};

% Select a subset of graph theory metrics to analyze.
allGTs = fieldnames(gt);
gtidxall = [7, 6, 5, 4, 15, 9, 13, 11];
allGTs = allGTs(gtidxall);

% Prepare time labels for the plots.
labelList = cell(size(stages_ptz, 1), 1);
for i = 1:size(stages_ptz, 1)
    val1 = stages_ptz(i, 1) / 2 - 300;
    val2 = stages_ptz(i, end) / 2 - 300;
    labelList{i} = mat2str(val1);
end
labelList{1} = 'Pre PTZ';
labelList{2} = 'Add PTZ';

% Initialize matrices to store statistical results.
meanMat = zeros(size(gt.gtom, 2), length(allGTs));
txtMat = cell(size(gt.gtom, 2), length(allGTs));
pMat = zeros(size(gt.gtom, 2), length(allGTs));

% 2. STATISTICAL ANALYSIS
for gtid = 1:length(allGTs)
    % Select and normalize the current graph metric.
    gtmat = gt.(allGTs{gtid});
    gtmat = (gtmat - min(gtmat(:))) / (max(gtmat(:)) - min(gtmat(:)));

    % Perform two-way repeated measures ANOVA for WT and MUT over time.
    [statoutput, AT] = ranova2(gtmat);
    pMat(:, gtid) = statoutput.pValue(2:2:end);

    % Compare genotypes for each stage.
    for jid = 1:size(gtmat, 2)
        data1 = gtmat(1:length(WT_fish), jid);
        data2 = gtmat(length(WT_fish) + 1:end, jid);

        % Calculate the normalized mean difference.
        meanMat(jid, gtid) = (nanmean(data1) - nanmean(data2)) / (nanmean(data1) + nanmean(data2));

        % Determine significance stars based on p-value.
        p = pMat(jid, gtid);
        sig_star = '';
        if p < 0.05 && p >= 0.01
            sig_star = '*';
        elseif p < 0.01 && p >= 0.001
            sig_star = '**';
        elseif p < 0.001
            sig_star = '***';
        end
        txtMat{jid, gtid} = sig_star;
    end
end

% 3. VISUALIZATION (BUBBLE PLOT)
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.4], ...
       'Color', 'w', 'Renderer', 'painters');
hold on;

% Define color map.
cols = redblueRealTecplot(256);
colormap(cols);
clim([-0.25, 0.25]);

% Plot the bubble chart.
for gtid = 1:length(allGTs)
    for i = [1, 4:size(gtmat, 2)]
        % Circle size represents statistical significance.
        circlesize = -log(pMat(i, gtid)) * 250;
        
        % Circle color represents the magnitude and direction of the difference.
        scatter(i * 2 - 1, gtid, circlesize, -meanMat(i, gtid), 'filled', ...
                'MarkerEdgeColor', 'none', 'LineWidth', 0.5);
        
        % Add black outline for significant points.
        if ~isempty(txtMat{i, gtid})
            scatter(i * 2 - 1, gtid, circlesize, [0, 0, 0], 'LineWidth', 2.5);
        end
    end
end

% Format the axes and add labels.
xlim([0, i * 2]);
ylim([0, gtid + 1]);
xticks(1:6:2 * size(gtmat, 2));
xticklabels(labelList(1:3:end));
yticks(1:length(allGTs));
yticklabels(allGTs);

xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Graph Theory Metric');
cb = colorbar;
cb.Label.FontSize = 12;
cb.Label.FontWeight = 'bold';
cb.Ticks = [-0.25, -0.15, 0, 0.15, 0.25];
cb.Label.String = 'Normalized Mean Difference (WT - MUT)';

%% Panel 2
%
% This script visualizes and statistically compares various graph theory
% metrics between WT and mutant fish across a 24-stage time course.
% It plots individual data points, mean trends, and annotates each plot
% with key statistical measures (Cohen's f and LME p-value).

%% 1. INITIALIZATION AND DATA PREPARATION

% Pre-process the 'gtom' data from a cell array into a matrix.
gt1 = gt;
gt1.gtom = zeros(35, 24);
for i = 1:35
    for j = 1:24
        gt1.gtom(i, j) = nanmean(gt.gtommap{i, j}(:, 4));
    end
end

% Reorder the fields of the 'gt' structure for consistent plotting.
orderedFields = fieldnames(gt1);
orderedFields = orderedFields([12, 7, 5, 15, 13, 6, 4, 9, 11, 14, 2, 1, 3, 8, 10]);
gt = orderfields(gt1, orderedFields);

% Get the relevant graph theory field names, excluding maps and other non-data fields.
rawFields = fieldnames(gt);
fieldNames = rawFields(cellfun(@(f) ~contains(f, 'map'), rawFields));
fieldNames = fieldNames(2:end-2);
nFields = numel(fieldNames);

% 2. PLOTTING AND STATISTICAL ANALYSIS
figure('Position', [0.1490, 0.1897, 1.7180, 0.5740] * 1000, 'Renderer', 'painters');

% Loop through each graph theory metric to create a subplot.
for i = 1:nFields
    subplot(2, 4, 9 - i); % Create a 2x4 subplot grid.
    hold on;
    
    % Set log scale for specific plots (if needed).
    if i == 8 || i == 7
        set(gca, 'YScale', 'log');
    end
    
    fieldData = gt.(fieldNames{i});
    data = fieldData;
    
    % Plot individual fish data points.
    for kwid = 1:size(data, 2)
        if ~any([2, 3] == kwid)
            data1 = data(1:length(WT_fish), kwid);
            data2 = data(length(WT_fish) + 1:end, kwid);
        else
            data1 = nan;
            data2 = nan;
        end
        h1 = scatter(kwid, data1, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.1, 'DisplayName', 'fish (WT)');
        h2 = scatter(kwid, data2, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.1, 'DisplayName', 'fish (MUT)');
        if kwid ~= 1
            set(h1, 'HandleVisibility', 'off');
            set(h2, 'HandleVisibility', 'off');
        end
    end

    % Calculate and plot mean trends and standard error of the mean (SEM).
    plotdata1 = nanmean(data(1:length(WT_fish), :));
    plotdata2 = nanmean(data(length(WT_fish) + 1:end, :));
    SEM1 = nanstd(data(1:length(WT_fish), :)) / sqrt(size(data, 1));
    SEM2 = nanstd(data(length(WT_fish) + 1:end, :)) / sqrt(size(data, 1));

    x = 1:size(data, 2);
    plot(x, plotdata1, 'b*-', 'LineWidth', 2);
    plot(x, plotdata2, 'r*-', 'LineWidth', 2);
    
    % Fill shaded areas for SEM.
    fill([x, fliplr(x)], [plotdata1 + SEM1, fliplr(plotdata1 - SEM1)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    fill([x, fliplr(x)], [plotdata2 + SEM2, fliplr(plotdata2 - SEM2)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % Perform statistical tests using Linear Mixed-Effects Model (LME).
    nSubjects = size(data, 1);
    nTimePoints = size(data, 2);
    Subject = repelem((1:nSubjects)', nTimePoints);
    Time = repmat((1:nTimePoints)', nSubjects, 1);
    Y = data';
    Y = Y(:);
    tbl = table(Subject, Time, Y);
    lme = fitlme(tbl, 'Y ~ Time + (Time|Subject)');
    anovaTable = anova(lme);
    pValue = anovaTable.pValue(2); % p-value for the 'Time' fixed effect.
    
    % Calculate Cohen's f for effect size using repeated measures ANOVA.
    T_rm = array2table(fieldData, 'VariableNames', strcat("Cond", string(1:nTimePoints)));
    T_rm.Fish = (1:size(fieldData, 1))';
    rm = fitrm(T_rm, sprintf('%s-%s ~ 1', T_rm.Properties.VariableNames{1}, T_rm.Properties.VariableNames{end}), ...
              'WithinDesign', table(strcat('t', string(1:nTimePoints))', 'VariableNames', {'Time'}));
    ranovaOutput = ranova(rm);
    SS_effect = ranovaOutput.SumSq(1);
    SS_total = sum(ranovaOutput.SumSq);
    eta2 = SS_effect / SS_total;
    cohen_f = sqrt(eta2 / (1 - eta2));
    
    % Add shaded regions for significant differences (p < 0.05).
    statoutput = ranova2(data);
    skipwin = 1;
    for kwid = 1:skipwin:size(data, 2) - skipwin
        p = statoutput.pValue(kwid * 2);
        if p < 0.05 && ~any([2, 3] == kwid)
            yl = ylim;
            p_patch = patch([kwid, kwid, kwid + skipwin, kwid + skipwin] - 0.5, [yl(2), yl(1), yl(1), yl(2)], 'k');
            if p < 0.01
                p_patch.FaceAlpha = 0.4;
            else
                p_patch.FaceAlpha = 0.15;
            end
            p_patch.EdgeColor = 'none';
        end
    end
    
    % Annotate the plot with statistical metrics.
    text(0.95, 0.95, sprintf('Cohen''s f = %.2g; LME(t) p = %.2g', cohen_f, pValue), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.2, 0.2, 0.6]);

    % Format the subplot.
    xlim([0, size(data, 2)]);
    ytickformat('%.2g');
    ylabel(fieldNames{i}, 'Interpreter', 'none');
    
    % Create time labels.
    timelabels = cell(1, size(data, 2));
    winsize = 300;
    slidewin = winsize / 2;
    for kwid = 1:size(data, 2)
        starttime = (kwid - 1) * slidewin;
        timelabels{kwid} = mat2str(starttime / 2);
    end
    timelabels{1} = 'Pre PTZ';
    timelabels{2} = 'Add PTZ';
    xticks(1:3:size(data, 2));
    xticklabels(timelabels(1:3:end));
end

% Add a common xlabel for the entire figure.
xlabel('Time (s)');

%% Panel C
% This script analyzes brain connectivity data from zebrafish to compare
% a wild-type group (WT) and a mutant group (SCN1).
% It performs a statistical analysis on node degree values in different brain
% regions and generates a significance plot.

% 1. Setup and Data Loading
clear all;
clc;
close all;

% Define fish groups
WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};
all_fish = [WT_fish, SCN1_fish];

% Load brain masks and ground truth data
run("zbrainmaskupdate.m");
Zbrain_Masks = Zbrain_shortlist;
load('whole_GT24.mat', 'gt');

% Define ground truth ID and stages
gt_id = 3;
all_gt_fields = fieldnames(gt);
selected_gt = all_gt_fields{gt_id};

% 2. Process Brain Regions
% Call the function to perform the statistical analysis.
% This separates the core logic from the main script.
[p_mat, mean_mat, subregion_names] = process_brain_regions(Zbrain_Masks, all_fish, gt, selected_gt);

% 3. Plotting Results
% Define plotting parameters
subregion_idx = find(~cellfun(@isempty, subregion_names)); % Find non-empty entries
subregions_to_plot = Zbrain_Masks(subregion_idx, 2);

% Call the function to generate the plot
plot_brain_significance(mean_mat, p_mat, subregions_to_plot, Zbrain_Masks, Zbrain_shortlist);

function [p_mat_all, mean_mat_all, subregion_names] = process_brain_regions(zbrain_masks, all_fish, gt_data, selected_gt)
% PROCESS_BRAIN_REGIONS
% Analyzes node degree for each brain subregion across different fish and stages.
%
% Inputs:
%   zbrain_masks: A cell array containing brain region information.
%   all_fish: A cell array of fish IDs.
%   gt_data: Ground truth data structure.
%   selected_gt: The name of the ground truth field to analyze.
%
% Outputs:
%   p_mat_all: A matrix of p-values from statistical tests.
%   mean_mat_all: A matrix of mean differences between groups.
%   subregion_names: A cell array of subregion names.

    % Initialize output variables
    num_regions = size(zbrain_masks, 1);
    num_stages = 24;
    p_mat_all = nan(num_regions, num_stages);
    mean_mat_all = zeros(num_regions, num_stages);
    subregion_names = cell(num_regions, 1);
    
    % Get all brain regions from a template
    per_region = getPerBrainRegionsAll(zbrain_masks, [0 0 0]);
    all_regions = fieldnames(per_region);
    
    count = 1;
    for regi = 1:length(all_regions)
        all_sub_regions = fieldnames(per_region.(all_regions{regi}));
        for subi = 1:length(all_sub_regions)
            
            % Skip if subregion has too few nodes
            sub_idx = per_region.(all_regions{regi}).(all_sub_regions{subi}).idx;
            if length(sub_idx) <= 10
                fprintf('Skipping region with <= 10 nodes: %s:%s\n', ...
                    all_regions{regi}, all_sub_regions{subi});
                continue;
            end
            
            subregion_names{count} = all_sub_regions{subi};
            
            for kwid = 1:num_stages
                data = [];
                geno_idx = [];
                fish_idx = [];
                
                % Loop through each fish to collect data for the current region and stage
                for fi = 1:length(all_fish)
                    % Load N (node count) for each fish
                    load(['reginfo/regioninfo_' all_fish{fi} '.mat'], 'N');
                    
                    % Calculate average node degree for the subregion
                    node_degree_avg = gt_data.(selected_gt){fi, kwid}(sub_idx, 4) / (N(fi) - 1);
                    
                    % Collect data for statistical test
                    data = [data; node_degree_avg];
                    
                    % Assign genotype and fish ID
                    if fi <= 17 % Assuming first 17 fish are WT
                        geno_idx = [geno_idx; ones(size(node_degree_avg))];
                    else % The rest are SCN1
                        geno_idx = [geno_idx; ones(size(node_degree_avg)) * 2];
                    end
                    fish_idx = [fish_idx; ones(size(node_degree_avg)) * fi];
                end
                
                % Perform statistical analysis (Linear Mixed-Effects model)
                p_mat_all(count, kwid) = getSignificant(data, geno_idx, fish_idx);
                
                % Calculate and store mean difference
                mean_wt = mean(data(geno_idx == 1));
                mean_scn1 = mean(data(geno_idx == 2));
                mean_mat_all(count, kwid) = (mean_wt - mean_scn1) / mean_wt;
            end
            count = count + 1;
        end
    end
end

function plot_brain_significance(mean_mat, p_mat, subregion_names, Zbrain_Masks, Zbrain_shortlist)
% PLOT_BRAIN_SIGNIFICANCE
% Generates a circle plot visualizing statistical significance and mean differences.
%
% Inputs:
%   mean_mat: Matrix of mean differences to be visualized.
%   p_mat: Matrix of p-values to determine circle size.
%   subregion_names: Cell array of subregion names for y-axis labels.
%   Zbrain_Masks: Full brain mask data for region colors.
%   Zbrain_shortlist: Another brain mask for group lines.

    figure('Renderer','painters');
    hold on;
    
    % Setup colors for regions
    reg_colors = getColours;
    reg_colors = flip(reg_colors);
    whole_idx = find(ismember(Zbrain_shortlist(:,2),' '));
    
    % Prepare data for plotting
    num_subregions = size(mean_mat, 1);
    num_stages = size(mean_mat, 2);
    
    % Create a color map based on sorted mean values
    even_size_mean_mat = -[mean_mat -mean_mat];
    [~, cols_idx] = sort(even_size_mean_mat(:));
    
    cols = redblueRealTecplot(length(cols_idx));
    col_mat = reshape(cols_idx, size(even_size_mean_mat));
    col_mat = col_mat(1:num_subregions, :);
    
    y_label_loc = [];
    count = 0;
    prev_region_group_idx = 1;
    
    for jid = 1:num_subregions
        % Determine region group for coloring and separators
        region_group_idx = sum(subregion_names{jid} < whole_idx);
        
        for kwid = 1:num_stages
            circle_size = -log(p_mat(jid, kwid)) * 20;
            circle_size = min(circle_size, 300); % Cap the circle size
            
            % Scatter plot for significance
            scatter(kwid*2-1, num_subregions - count, circle_size, ...
                cols(col_mat(jid, kwid), :), 'filled');
            
            % Add a black border for significant results
            if p_mat(jid, kwid) < 0.05
                scatter(kwid*2-1, num_subregions - count, circle_size, [0 0 0], 'LineWidth', 2);
            end
        end
        
        % Add region name label to the plot
        if ~isinf(circle_size)
            count = count + 1;
            y_label_loc = [y_label_loc; count];
            region_name = subregion_names{jid};
            region_name = strrep(region_name, ' enriched area', '');
            region_name = strrep(region_name, ' Enriched Area', '');
            text(-0.2, num_subregions - count + 1, region_name, ...
                'Color', reg_colors{region_group_idx}, 'HorizontalAlignment', 'right');
        end
        
        % Add horizontal line to separate major brain regions
        if prev_region_group_idx ~= region_group_idx && region_group_idx < 4
            yline(num_subregions - count + 1.5, '--');
            prev_region_group_idx = region_group_idx;
        end
    end
    
    % Final plot adjustments
    ylim([0 num_subregions + 1]);
    xticks(1:2:num_stages*2);
    xlabel('Time (s)');
    ylabel('Brain Subregion');
    title('Genotype Difference in Node Degree');
    yticklabels([]);
end