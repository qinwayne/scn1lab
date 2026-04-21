%% Panel A

%% 1. INITIALIZATION AND DATA PREPARATION
clear all;
clc;
close all;

% Load the data file.
load('whole_GT24.mat');

winsz = 300;
slidsz = winsz/2;
stages_ptz(1,:) = 100:100+winsz;

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
    try
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

%% Panel C
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
num_stages = 24;

% We calculate and store the region maps for all fish upfront so the loop runs fast
fprintf('Precomputing brain regions for all fish...\n');
all_fish_regions = cell(length(all_fish), 1);

for fi = 1:length(all_fish)
    % Extract ROI locations for the current fish
    % Note: using (:, 1:3) assuming data is N rows (neurons) and 3 columns (X,Y,Z).
    % If your data is 3 rows and N columns, change this to gt.degreemap{fi, 1}(1:3, :)'
    fish_ROIs = gt.degreemap{fi, 1}(:, 1:3); 
    
    % Map this specific fish's ROIs to the Zbrain atlas
    all_fish_regions{fi} = getPerBrainRegionsAll(Zbrain_Masks, fish_ROIs);
end
fprintf('Precomputation complete.\n');

% Use Fish 1 as the structural reference to get region names
all_regions = fieldnames(all_fish_regions{1});

% Preallocate arrays for speed 
max_regions = 200; 
p_mat = nan(max_regions, num_stages);
mean_mat = zeros(max_regions, num_stages);
subregion_names = cell(max_regions, 1);

% 2. Data Extraction and Statistical Analysis Loop
count = 1;

for regi = 1:length(all_regions)
    all_sub_regions = fieldnames(all_fish_regions{1}.(all_regions{regi}));
    
    for subi = 1:length(all_sub_regions)
        
        % Check if this region is viable by looking at the first fish
        % (Skip if the subregion has too few nodes to be statistically relevant)
        ref_sub_idx = all_fish_regions{1}.(all_regions{regi}).(all_sub_regions{subi}).idx;
        if length(ref_sub_idx) <= 10
            continue;
        end
        
        subregion_names{count} = all_sub_regions{subi};
        
        for kwid = 1:num_stages
            data = [];
            geno_idx = [];
            fish_idx = [];
            
            % Loop through each fish
            for fi = 1:length(all_fish)
                
                % ---> Extract the UNIQUE sub_idx for THIS specific fish <---
                fish_sub_idx = all_fish_regions{fi}.(all_regions{regi}).(all_sub_regions{subi}).idx;
                
                % Only process if this fish actually has neurons in this region
                if isempty(fish_sub_idx)
                    continue;
                end
                
                % Load node count (N) for current fish
                load(['reginfo/regioninfo_' all_fish{fi} '.mat'], 'N');
                
                % Calculate node degree using this fish's specific sub_idx
                node_degree_avg = gt.(selected_gt){fi, kwid}(fish_sub_idx, 4) / (N(fi) - 1);
                data = [data; node_degree_avg];
                
                % Assign genotype: 1 for WT (first 17 fish), 2 for SCN1 (mutant)
                if fi <= 17 
                    geno_idx = [geno_idx; ones(size(node_degree_avg))];
                else 
                    geno_idx = [geno_idx; ones(size(node_degree_avg)) * 2];
                end
                fish_idx = [fish_idx; ones(size(node_degree_avg)) * fi];
            end
            
            % Stats Calculation (Ensure there is data before testing to prevent crashes)
            if ~isempty(data)
                p_mat(count, kwid) = getSignificant(data, geno_idx, fish_idx);
                
                % Mean Difference Calculation (Fixes the NaN propagation bug)
                mean_wt = mean(data(geno_idx == 1), 'omitnan');
                mean_scn1 = mean(data(geno_idx == 2), 'omitnan');
                mean_mat(count, kwid) = (mean_wt - mean_scn1) / mean_wt;
            end
        end
        count = count + 1;
    end
end

% Truncate matrices to remove trailing NaNs/Zeros
num_valid_regions = count - 1;
p_mat = p_mat(1:num_valid_regions, :);
mean_mat = mean_mat(1:num_valid_regions, :);
subregion_names = subregion_names(1:num_valid_regions);

% 3. Plotting the Bubble Chart
figure('Renderer', 'painters');
hold on;

% Setup colors and indices
reg_colors = flip(getColours());
whole_idx = find(ismember(Zbrain_shortlist(:, 2), ' '));

% Create a colormap based on sorted mean values
even_size_mean_mat = -[mean_mat, -mean_mat];
[~, cols_idx] = sort(even_size_mean_mat(:));
cols = redblueRealTecplot(length(cols_idx));

col_mat = reshape(cols_idx, size(even_size_mean_mat));
col_mat = col_mat(1:num_valid_regions, :);

plot_count = 0;
prev_region_group_idx = 1;

for jid = 1:num_valid_regions
    % Determine region group for text coloring and separators
    region_group_idx = sum(subregion_names{jid} < whole_idx);
    
    for kwid = 1:num_stages
        % Calculate circle size based on p-value
        circle_size = -log(p_mat(jid, kwid)) * 20;
        circle_size = min(circle_size, 300); % Cap size to prevent massive overlap
        
        % Scatter plot for the bubble
        scatter(kwid * 2 - 1, num_valid_regions - plot_count, circle_size, ...
            cols(col_mat(jid, kwid), :), 'filled');
        
        % Add black outline if statistically significant (p < 0.05)
        if p_mat(jid, kwid) < 0.05
            scatter(kwid * 2 - 1, num_valid_regions - plot_count, circle_size, [0 0 0], 'LineWidth', 2);
        end
    end
    
    % Add text label to Y-Axis
    if ~isinf(circle_size) && ~isnan(circle_size)
        plot_count = plot_count + 1;
        
        % Clean up the region name string
        region_name = subregion_names{jid};
        region_name = strrep(region_name, ' enriched area', '');
        region_name = strrep(region_name, ' Enriched Area', '');
        
        text(-0.2, num_valid_regions - plot_count + 1, region_name, ...
            'Color', reg_colors{region_group_idx}, 'HorizontalAlignment', 'right');
    end
    
    % Draw horizontal separator line between major brain groups
    if prev_region_group_idx ~= region_group_idx && region_group_idx < 4
        yline(num_valid_regions - plot_count + 1.5, '--');
        prev_region_group_idx = region_group_idx;
    end
end

% Final Plot Adjustments
ylim([0, num_valid_regions + 1]);
xticks(1:2:num_stages * 2);
xlabel('Time (s)');
ylabel('Brain Subregion');
title('Genotype Difference in Node Degree');
yticklabels([]); % Hide standard y-ticks since we use text()