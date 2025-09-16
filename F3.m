%% Panel B
% This script compares seizure event density maps between WT and SCN1
% mutant fish. It performs a voxel-wise statistical analysis and visualizes
% the results as a series of 2D slices.

% 1. INITIALIZATION & DATA LOADING
clear all;
load("density_comparison2.mat");

% Define parameters for the analysis.
pixel_add = 1; % Defines the size of the local neighborhood (e.g., a 3x3x3 cube) for the statistical test.
finalMat = zeros(targetSize); % Initialize a 3D matrix to store the results.
all_fish = [WT_fish, SCN1_fish]; % Combine fish lists for easier iteration.

% 2. STATISTICAL COMPARISON (VOXEL-WISE)
% This nested loop structure iterates through each voxel in the 3D space,
% creating a local neighborhood and performing a statistical test.
for xi = 1 + pixel_add : targetSize(1) - pixel_add
    for yi = 1 + pixel_add : targetSize(2) - pixel_add
        for zi = 1 + pixel_add : targetSize(3) - pixel_add
            % Extract voxel data from the local neighborhood for WT and MUT fish.
            wt_mov_den = [];
            mut_mov_den = [];

            for fi = 1 : length(WT_fish)
                % Extract a local cube of data for the current WT fish.
                data = reszden{fi}(xi - pixel_add : xi + pixel_add, ...
                                   yi - pixel_add : yi + pixel_add, ...
                                   zi - pixel_add : zi + pixel_add);
                wt_mov_den = [wt_mov_den; data(:)]; % Reshape and append.
            end

            for fi = length(WT_fish) + 1 : length(all_fish)
                % Extract a local cube of data for the current MUT fish.
                data = reszden{fi}(xi - pixel_add : xi + pixel_add, ...
                                   yi - pixel_add : yi + pixel_add, ...
                                   zi - pixel_add : zi + pixel_add);
                mut_mov_den = [mut_mov_den; data(:)]; % Reshape and append.
            end
            
            % Handle cases where one or both groups have no data.
            % The `try...end` block prevents the code from crashing.
            try
                % Perform the Mann-Whitney U test (ranksum) to compare the
                % two groups' distributions. This is a non-parametric test.
                [p, h] = ranksum(wt_mov_den, mut_mov_den);

                % If the p-value is significant (p < 0.05), store the
                % difference in mean density in the `finalMat`.
                if p < 0.05
                    finalMat(xi, yi, zi) = nanmean(wt_mov_den) - nanmean(mut_mov_den);
                end
            end
        end
    end
end

% 3. VISUALIZATION OF RESULTS
% Generate a multi-panel figure to display the significant differences
% across different slices of the brain.
close all;
fig = figure('position', 1000 * [1.0000, 1.0490, 0.5600, 0.1887], 'Renderer', 'painters');

% Apply a custom colormap.
colormap(flip(redblueRealTecplot));

% Define scaling factors for overlaying a reference image (e.g., a brain atlas).
scalex = (850 - 0) / targetSize(1);
scaley = (430 - 0) / targetSize(2);

count = 1;
% Loop through the z-axis (slices) of the `finalMat`.
for i = 1 : 1 : size(finalMat, 3) - 1
    subplot(4, 4, count);
    
    % Prepare the data for display.
    data_temp = squeeze(rot90(finalMat(:, :, i), 2));
    
    % Create an alpha mask to make non-significant voxels (value 0) transparent.
    imAlpha = ones(size(data_temp));
    imAlpha(data_temp == 0) = 0;
    
    % Use `imagesc` to display the data with the alpha mask.
    imagesc(data_temp, 'AlphaData', imAlpha);
    
    % Set the color axis limits to a fixed range for consistent visualization.
    clim([-0.3, 0.3]);
    
    % Format the axes to be clean and minimal.
    set(gca, 'box', 'off', 'XTickLabel', [], 'XTick', [], 'YTickLabel', [], ...
             'YTick', [], 'TickDir', 'out', 'xcolor', [1, 1, 1], ...
             'ycolor', [1, 1, 1]);

    % Overlay a reference line (e.g., a midline or a specific structure).
    hold on;
    plot(tx / scaley - 20, 220 - ty / scalex, 'k-.');
    
    count = count + 1;
end

% 4. FINAL FIGURE FORMATTING
% Create a new axes object for a common colorbar and labels.
h = axes(fig, 'visible', 'off');
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';

% Add a color bar.
handleToColorBar = colorbar;
set(handleToColorBar, 'YTickLabel', []);
handleToColorBar.Location = 'west';

% Add custom labels to the colorbar.
hYLabel = ylabel(handleToColorBar, {'MUT', 'WT'});
handleToColorBar.Position = [0.92, 0.1349, 0.0381, 0.7651];
%% Panel C
% This script generates a figure with multiple subplots to visualize the
% maximum proportion of synchronized ROIs per brain region for WT and
% SCN1 mutant fish.

% 1. INITIALIZE FIGURE
% Create a new figure with a specified size.
fig = figure('Position', [1000, 944.3, 560, 293.3]);

% 2. LOOP THROUGH BRAIN REGIONS AND PLOT DATA
% The loop iterates through four brain regions (i = 1 to 4).
for i = 1:4
    % Use subplots to arrange the plots.
    % The first two regions are plotted in the top two panels, and the
    % last two regions are plotted in the bottom two panels.
    if i > 2
        subplot(4, 1, [1:2]);
    else
        subplot(4, 1, [3:4]);
        % Add a legend only to the bottom subplot to avoid redundancy.
        legend('Location', 'best');
    end
    
    hold on;
    xlim([0, 11]);
    
    % Calculate the maximum ratio of synchronized ROIs for each fish.
    % The ratio is the maximum number of synchronized ROIs in a given
    % region divided by the total number of ROIs for that fish.
    WTratio = zeros(length(WT_fish), 4);
    MUTratio = zeros(length(SCN1_fish), 4);
    
    for fi = 1:length(all_fish)
        if fi <= length(WT_fish)
            % Calculate WT ratio for the current region.
            WTratio(fi, i) = max(WTidx(i, :, fi), [], 2) / N{fi};
        else
            % Calculate mutant (MUT) ratio for the current region.
            % Note the index for MUTidx and N is adjusted.
            mutant_idx = fi - length(WT_fish);
            MUTratio(mutant_idx, i) = max(MUTidx(i, :, mutant_idx), [], 2) / N{fi};
        end
    end
    
    % Plot individual data points as scatter plots.
    % `HandleVisibility`, when set to 'off', prevents these points from appearing in the legend.
    scatter((i*2+1) - 1, WTratio(:, i), 'Marker', '.', 'MarkerEdgeColor', 'blue', 'HandleVisibility', 'off');
    scatter((i*2+1), MUTratio(:, i), 'Marker', '.', 'MarkerEdgeColor', 'red', 'HandleVisibility', 'off');
    
    % Create boxcharts to visualize the data distribution.
    if i == 1
        % For the first region (i=1), plot boxcharts with legend entries.
        boxchart(i*2 + ones(length(WT_fish), 1) - 1, WTratio(:, i), 'BoxFaceColor', 'b', 'MarkerStyle', '.', 'MarkerColor', 'b', 'DisplayName', 'WT');
        boxchart(i*2 + ones(length(SCN1_fish), 1), MUTratio(:, i), 'BoxFaceColor', 'r', 'MarkerStyle', '.', 'MarkerColor', 'r', 'DisplayName', 'MUT');
    else
        % For subsequent regions, plot boxcharts without legend entries.
        boxchart(i*2 + ones(length(WT_fish), 1) - 1, WTratio(:, i), 'BoxFaceColor', 'b', 'MarkerStyle', '.', 'MarkerColor', 'b', 'HandleVisibility', 'off');
        boxchart(i*2 + ones(length(SCN1_fish), 1), MUTratio(:, i), 'BoxFaceColor', 'r', 'MarkerStyle', '.', 'MarkerColor', 'r', 'HandleVisibility', 'off');
    end
    
    % Plot a line connecting the mean values of WT and MUT for each region.
    plot([i*2, i*2+1], [nanmean(WTratio(:, i)), nanmean(MUTratio(:, i))], 'k', 'HandleVisibility', 'off');
    
    % Add significance stars based on a previous statistical test.
    % The `h` variable is a boolean indicating significance, and `sig_star` holds the annotation string.
    if ~isnan(h) && h
        plot([i*2, i*2+1], [0.5, 0.5], 'k', 'HandleVisibility', 'off');
        text(i*2 + 0.4, 0.52, sig_star, 'FontSize', 15);
    end
    
    % Turn off the X-axis for the current subplot.
    current_axis = gca;
    current_axis.XAxis.Visible = 'off';
end

% 3. FINAL PLOT FORMATTING
% Re-enable and format the X-axis on the bottom subplot.
subplot(4, 1, [3:4]);
current_axis = gca;
current_axis.XAxis.Visible = 'on';

% Define X-axis ticks and labels.
xticks(2.5:2:9);
xticklabels([]); % Initially hide tick labels.

% Use a custom function to get region names and colors.
[color_labels, RegionList] = getColours();
% The `getColours` function is a placeholder and should be provided separately.

% Add text labels for each region at the bottom of the plot.
for i = 1:4
    % Add the first four letters of each region name as a label.
    text(2*i + 1.1, -0.012, RegionList{i}(1:4), 'Color', color_labels{i}, 'rotation', 0, 'HorizontalAlignment', 'right', 'FontSize', 13);
end

% Create a common Y-label for the entire figure.
han = axes(fig, 'visible', 'off');
han.YLabel.Visible = 'on';
ylabel(han, 'Proportion of total ROIs');
%% Panel D
% This script identifies and plots subregions where the proportion of
% synchronized ROIs significantly differs between WT and SCN1 mutant fish.
% It calculates this proportion by normalizing the number of active ROIs
% by the total number of ROIs for each animal.

% 1. INITIALIZATION & DATA LOADING
% Clear previous plots and load data.
close all;
load("stat_number_region_shortlist.mat");

% Initialize counters and figure.
count = 0; % Counter for the number of significant plots.
fig = figure;
hold on;

% Get region colors and names from a custom function.
[collabels, RegionList] = getColours();

% 2. PLOTTING
% Loop through each brain region and subregion to perform statistical tests.
for regi = 1:size(WTidx, 1)
    for subi = 1:size(WTidx, 2)
        % Calculate the proportion of synchronized ROIs for each fish.
        % The number of ROIs is normalized by the total number of ROIs (N)
        % for each animal to account for variations in ROI count.
        nowt = squeeze(WTidx(regi, subi, :))' ./ N(1:length(WT_fish));
        nomut = squeeze(MUTidx(regi, subi, :))' ./ N(1+length(WT_fish):end);

        % Calculate the mean change percentage for reference (not used for plotting).
        change_perc = (nanmean(nomut) - nanmean(nowt)) / nanmean(nowt) * 100;

        % Only proceed if both WT and mutant groups have non-zero data.
        if sum(nowt == 0) == 0 && sum(nomut == 0) == 0

            % If the p-value is significant (p < 0.05), plot the results.
            if p(regi,subi) < 0.05
                count = count + 1;

                % Plot individual data points for the mutant group.
                % The data is normalized by the mean of the WT group and
                % multiplied by 100 to represent a percentage of the WT mean.
                scatter(ones(length(SCN1_fish), 1) * count, nomut / nanmean(nowt) * 100, ...
                        'Marker', '.', 'MarkerEdgeColor', collabels{regi}, 'HandleVisibility', 'off');

                % Plot boxchart for the mutant group.
                h_mut = boxchart(ones(length(SCN1_fish), 1) * count, nomut / nanmean(nowt) * 100, ...
                                 'BoxFaceColor', collabels{regi}, 'MarkerStyle', '.', ...
                                 'MarkerColor', collabels{regi}, 'DisplayName', 'MUT', 'HandleVisibility', 'off');

                % Draw a horizontal line at 100% to represent the WT mean.
                yline(100, '-.');

                % Set the handle visibility for the first significant plot
                % to 'on' so it appears in the legend.
                if count == 1
                    set(h_mut, 'HandleVisibility', 'on');
                end
                
                % Store the name of the significant subregion.
                temp_name = MUTnames{regi, subi};
                temp_name = strrep(temp_name, '_', ' ');
                temp_name = strrep(temp_name, 'cephalon', '.');
                xlabelnames{count} = temp_name;
                
                % Store the indices for later use (e.g., in other scripts).
                subregidx(count, :) = [regi, subi];
            end
        end
    end

    count = count + 1;
    WTratio = max(squeeze(WTidx(regi, :, :))) ./ N(1:length(WT_fish))';
    MUTratio = max(squeeze(MUTidx(regi, :, :))) ./ N(1+length(WT_fish):end)';

    % Plot scatter and boxchart for the main region.
    scatter(ones(length(SCN1_fish), 1) * count, MUTratio / nanmean(WTratio) * 100, ...
            'Marker', '.', 'MarkerEdgeColor', collabels{regi}, 'HandleVisibility', 'off');
    boxchart(ones(length(SCN1_fish), 1) * count, MUTratio / nanmean(WTratio) * 100, ...
             'BoxFaceColor', collabels{regi}, 'MarkerStyle', '.', ...
             'MarkerColor', collabels{regi}, 'HandleVisibility', 'off');
    xlabelnames{count} = RegionList{regi};

end

% 3. PLOT FINALIZATION
% Set X-axis tick locations, labels, and plot range.
xticks(1 : (count));
xlim([0.5, count + 0.5]);
xticklabels(xlabelnames);
xtickangle(45); % Consider rotating labels if they overlap.

%% Panel E
%
% This script computes and plots the histogram of pairwise correlations
% between ROIs for WT and SCN1 mutant fish. It compares these distributions
% to a randomized null model and across different stages of a PTZ experiment.

% 1. INITIALIZATION AND DATA PROCESSING
% Clear all variables and close figures from previous runs.
clear all;
close all;

% Initialize matrices to store correlation vectors for each group and condition.
AllcorrVecWT = [];
AllcorrVecMUT = [];
AllcorrVecRND = []; % For the null model.
skipsz = 10; % Downsampling factor for the time series data to speed up processing.

% Loop through each fish to calculate correlation matrices.
for fi = 1:length(all_fish)
    % Determine the genotype and file path based on the fish index.
    if fi > length(WT_fish)
        genotype = 'Hom';
        color = 'red';
    else
        genotype = 'WT';
        color = 'blue';
    end
    
    % Display current fish ID for tracking progress.
    fish_id = all_fish{fi};
    fprintf('Processing fish: %s (%s)\n', fish_id, genotype);
    
    % Construct the file path and load the data.
    filepath = fullfile('\The University of Melbourne\Research\SCN1\data\SCN1LabData', ...
                        genotype, ['raw_fish_std_fmt_' fish_id '.mat']);
    load(filepath, 'fish_stim_trains', 'allcond', 'pre_ptz', 'early_ptz', 'late_ptz');

    % Downsample the time series data.
    df_all = fish_stim_trains{1, 1}(1:skipsz:end, :);
    
    % --- Null Model: Randomized Data ---
    % Create a random permutation of the data to generate a null distribution.
    % This breaks any existing temporal correlations.
    SMx = reshape(randperm(size(df_all, 1) * size(df_all, 2)), size(df_all, 1), size(df_all, 2));
    df_rand = df_all(SMx);
    
    % Calculate the correlation matrix for the randomized data.
    corrMatNull = corrcoef(df_rand');
    corrMatNull(eye(size(corrMatNull)) == 1) = NaN; % Set diagonal (self-correlation) to NaN.
    AllcorrVecRND = [AllcorrVecRND, reshape(corrMatNull, 1, [])];
    
    % --- Experimental Conditions ---
    condcorrVec = [];
    % Loop through each condition (pre, early, late PTZ).
    for kwid = 1:3
        % Select the data corresponding to the current time window.
        keywords = allcond{kwid};
        time_indices = eval([keywords '_ptz']);
        df_kw = df_all(:, time_indices);
        
        % Calculate the correlation matrix.
        corrMat = corrcoef(df_kw');
        corrMat(eye(size(corrMat)) == 1) = NaN; % Set diagonal to NaN.
        
        % Store the reshaped correlation matrix for the current fish and condition.
        corrVec = reshape(corrMat, 1, []);
        condcorrVec = [condcorrVec; corrVec];
    end
    
    % Append the results to the appropriate group's matrix.
    if fi > length(WT_fish)
        AllcorrVecMUT = [AllcorrVecMUT, condcorrVec];
    else
        AllcorrVecWT = [AllcorrVecWT, condcorrVec];
    end
end

% 2. VISUALIZATION OF RESULTS
% Generate a figure with four subplots to display the histograms.
close all;
fig = figure('Position', [846, 151, 354, 714]);
linestyles = {'-', '-', '-'};
titlelists = {'Null Model', 'Pre-PTZ', 'Early PTZ: 100s-950s', 'Late PTZ: 950s-1800s'};

% --- Subplot 1: Null Model ---
subplot(4, 1, 1);
hold on;
% Plot histogram for the randomized data (null model).
histogram(AllcorrVecRND, 500, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', ...
          'FaceColor', 'none', 'LineStyle', linestyles{1}, 'Normalization', 'probability');
xlim([-0.3, 1]);
xticks([]);
yticks([]);
title(titlelists{1});

% --- Subplot 2: Pre-PTZ Condition ---
subplot(4, 1, 2);
hold on;
% Plot histograms for WT and MUT groups.
histogram(AllcorrVecWT(1, :), 500, 'DisplayStyle', 'stairs', 'EdgeColor', 'b', ...
          'FaceColor', 'none', 'LineStyle', linestyles{1}, 'Normalization', 'probability');
histogram(AllcorrVecMUT(1, :), 500, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', ...
          'FaceColor', 'none', 'LineStyle', linestyles{1}, 'Normalization', 'probability');
xlim([-0.3, 1]);
xticks([]);
yticks([]);
title(titlelists{2});

% --- Subplot 3: Early PTZ Condition ---
subplot(4, 1, 3);
hold on;
histogram(AllcorrVecWT(2, :), 500, 'DisplayStyle', 'stairs', 'EdgeColor', 'b', ...
          'FaceColor', 'none', 'LineStyle', linestyles{2}, 'Normalization', 'probability');
histogram(AllcorrVecMUT(2, :), 500, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', ...
          'FaceColor', 'none', 'LineStyle', linestyles{2}, 'Normalization', 'probability');
xlim([-0.3, 1]);
xticks([]);
yticks([]);
title(titlelists{3});

% --- Subplot 4: Late PTZ Condition ---
subplot(4, 1, 4);
hold on;
histogram(AllcorrVecWT(3, :), 500, 'DisplayStyle', 'stairs', 'EdgeColor', 'b', ...
          'FaceColor', 'none', 'LineStyle', linestyles{3}, 'Normalization', 'probability');
histogram(AllcorrVecMUT(3, :), 500, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', ...
          'FaceColor', 'none', 'LineStyle', linestyles{3}, 'Normalization', 'probability');
legend('WT', 'MUT', 'Location', 'best');
xlim([-0.3, 1]);
yticks([]);
xlabel('Correlation');
title(titlelists{4});

% --- Common Y-axis Label ---
% Create an invisible axes object to add a shared y-axis label for the entire figure.
han = axes(fig, 'visible', 'off');
han.YLabel.Visible = 'on';
ylabel(han, 'Probability Density Function');

%% Panel F

% This script visualizes and statistically analyzes the mean pairwise
% correlation of ROIs over time for WT and SCN1 mutant zebrafish.

% 1. INITIALIZATION & DATA LOADING
% Clear workspace and load pre-computed data.
close all;
clear all;
load("correlation.mat");
data = meancorrfishlist;

% Initialize figure for plotting.
figure('position', [1177, 811, 560, 377], 'Renderer', 'painters');
hold on;

% Initialize matrices to store mean data for plotting.
plotdata1 = zeros(1, size(data, 2));
plotdata2 = zeros(1, size(data, 2));

% 2. PLOT INDIVIDUAL FISH DATA AND CALCULATE MEANS
% Loop through each time window (column) of the data.
for kwid = 1:size(data, 2)
    data1 = data(1:length(WT_fish), kwid);
    data2 = data(length(WT_fish) + 1:end, kwid);

    % Plot individual data points as scatter plots.
    h1 = scatter(kwid, data1, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.1, 'DisplayName', 'fish (WT)');
    h2 = scatter(kwid, data2, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.1, 'DisplayName', 'fish (MUT)');
    
    % Only show the legend handle for the first scatter plot to avoid clutter.
    if kwid ~= 1
        set(h1, 'HandleVisibility', 'off');
        set(h2, 'HandleVisibility', 'off');
    end

    % Store the mean of each group for later plotting.
    plotdata1(:, kwid) = nanmean(data1);
    plotdata2(:, kwid) = nanmean(data2);
end

% Add the legend for the scatter plots.
legend('Location', 'best');

% 3. REPEATED MEASURES ANOVA
% This section performs a two-way repeated measures ANOVA.
clc;

% Create a table for the repeated measures model.
fishIDs = (1:size(meancorrfishlist, 1))';
condNames = strcat("Cond", string(1:size(meancorrfishlist, 2)));
T = array2table(meancorrfishlist, 'VariableNames', condNames);
T.Fish = fishIDs;

% Fit the repeated measures model.
rm = fitrm(T, sprintf('%s-%s ~ 1', condNames{1}, condNames{end}), ...
           'WithinDesign', table(condNames', 'VariableNames', {'Condition'}));

% Run the ANOVA.
ranovaOutput = ranova(rm);

% Calculate effect size (Cohen's f).
SS_effect = ranovaOutput.SumSq(1); % Sum of squares for the Condition effect.
SS_total = sum(ranovaOutput.SumSq);
eta2 = SS_effect / SS_total;
cohen_f = sqrt(eta2 / (1 - eta2));
fprintf('Cohen''s f: %.4f\n', cohen_f);

% 4. HIGHLIGHT SIGNIFICANT TIME POINTS
% This section iterates through the statistical output and highlights
% time points where there is a significant difference (p < 0.05).
skipwin = 1; % The step size for the loop.
count1 = 0; % Counter for p<0.05 legend.
count2 = 0; % Counter for p<0.01 legend.

for kwid = 1 : skipwin : size(data, 2) - skipwin
    % The `statoutput` variable contains p-values for each comparison.
    p = statoutput.pValue(kwid * 2);

    if p < 0.05
        % Create a shaded patch to highlight the significant time window.
        p1 = patch([kwid, kwid, kwid + skipwin, kwid + skipwin] - 0.5, ...
                   [0, 0.9, 0.9, 0], 'k');
        
        % Set the transparency and legend entry based on the p-value.
        if p < 0.05 && p > 0.01
            p1.FaceAlpha = 0.15;
            if count1 == 0
                set(p1, 'DisplayName', 'p<0.05');
                count1 = 1;
            else
                set(p1, 'HandleVisibility', 'off');
            end
        elseif p < 0.01
            p1.FaceAlpha = 0.4;
            if count2 == 0
                set(p1, 'DisplayName', 'p<0.01');
                count2 = 1;
            else
                set(p1, 'HandleVisibility', 'off');
            end
        end
        p1.EdgeColor = 'none';
    end
end

% 5. PLOT MEAN AND SEM OVER TIME
% Calculate the Standard Error of the Mean (SEM) for each group.
SEM1 = std(data(1:length(WT_fish), :), 'omitnan') / sqrt(length(WT_fish));
SEM2 = std(data(1 + length(WT_fish):end, :), 'omitnan') / sqrt(length(SCN1_fish));

% Use a moving mean to smooth the data for plotting.
movwin = 3;
x_smooth = movmean(1:size(data, 2), movwin, "omitnan");
y1_smooth = movmean(nanmean(plotdata1, 1), movwin, "omitnan");
y2_smooth = movmean(nanmean(plotdata2, 1), movwin, "omitnan");

% Plot the smoothed mean correlation lines.
plot(x_smooth, y1_smooth, 'b', 'LineWidth', 1, 'DisplayName', 'WT mean');
plot(x_smooth, y2_smooth, 'r', 'LineWidth', 1, 'DisplayName', 'MUT mean');

% Create shaded areas for the SEM.
x2 = [x_smooth, fliplr(x_smooth)];
% WT SEM
inBetween1 = [y1_smooth + SEM1, fliplr(y1_smooth - SEM1)];
fill(x2, inBetween1, 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
     'DisplayName', 'SEM (WT)');
% MUT SEM
inBetween2 = [y2_smooth + SEM2, fliplr(y2_smooth - SEM2)];
fill(x2, inBetween2, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
     'DisplayName', 'SEM (MUT)');

% Final plot formatting.
xlim([0, size(data, 2)]);
ylabel('Mean Correlation');
xlabel('Time (s)');

% Generate x-axis tick labels based on time windows.
kwid = 1;
while kwid * slidewin + winsize + 100 < 4200
    starttime = (kwid - 1) * slidewin + 100;
    endtime = kwid * slidewin + winsize + 100;
    if starttime < 600
        timelabels{kwid} = 'Pre PTZ';
    else
        timelabels{kwid} = [mat2str(starttime/2 - 300), '~', mat2str(endtime/2 - 300)];
    end
    kwid = kwid + 1;
end
xticks(3:10:size(data, 2));
xticklabels(timelabels(3:10:size(data, 2)));

% Update the legend to include all plotted elements.
legend('Location', 'best');