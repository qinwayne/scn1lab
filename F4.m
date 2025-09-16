%% Panel A, B

% This script visualizes the decay of pairwise neural correlation with
% increasing distance. It also compares the parameters of a power-law fit
% to these curves between WT and SCN1a mutant fish.

% 1. VISUALIZATION OF CORRELATION-DISTANCE RELATIONSHIP
close all;
clear all;
load('distance.mat');

% Create the main figure.
fig = figure;

% Define the x-axis for the plots (Euclidean distance).
x = (0:size(corrdisMat, 1) - 1) * distint * 3 / 2;

% Loop through the three experimental conditions.
for kwid = 1:3
    % Subplot for the overall correlation-distance curve.
    subplot(3, 2, 2 * kwid - 1);
    hold on;

    % Plot the WT mean correlation and SEM shaded area.
    y1 = nanmean(corrdisMat(:, 1:length(WT_fish), kwid), 2)';
    dy1 = nanstd(corrdisMat(:, 1:length(WT_fish), kwid), 0, 2)' / sqrt(length(WT_fish));
    curve1 = y1 + dy1;
    curve2 = y1 - dy1;
    inBetween = [curve1, fliplr(curve2)];
    plot(x, y1, 'b');
    fill([x, fliplr(x)], inBetween, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

    % Plot the MUT mean correlation and SEM shaded area.
    y2 = nanmean(corrdisMat(:, length(WT_fish) + 1:end, kwid), 2)';
    dy2 = nanstd(corrdisMat(:, length(WT_fish) + 1:end, kwid), 0, 2)' / sqrt(length(SCN1_fish));
    curve1 = y2 + dy2;
    curve2 = y2 - dy2;
    inBetween = [curve1, fliplr(curve2)];
    plot(x, y2, 'r');
    fill([x, fliplr(x)], inBetween, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

    % Set plot properties.
    ylim([0.03, 0.55]); % Adjusted based on condition.
    xlim([0, 900]);
    set(gca, 'YScale', 'log', 'color', 'none');

    % Add specific plot titles and formatting based on the condition.
    switch kwid
        case 1
            legend({'Mean (WT)', 'SEM (WT)', 'Mean (scn1lab^{-/-})', 'SEM (scn1lab^{-/-})'}, ...
                   'Box', 'off', 'position', [0.6263, 0.1070, 0.1601, 0.1401]);
            xticklabels([]);
            title('(i) Pre PTZ');
            ylim([0.03, 0.12]);
        case 2
            xticklabels([]);
            title('(ii) 50 ~ 950s');
            ylim([0.2, 0.4]);
        case 3
            title('(iii) 950 ~ 1800s');
            ylim([0.3, 0.55]);
    end

    % Add the enlarged subplots for Pre-PTZ condition.
    if kwid == 1
        % WT enlarged plot.
        subplot(3, 2, 2);
        hold on;
        plot(x, y1, 'b-.', 'LineWidth', 1);
        plot(x(7:13), y1(7:13), 'b-', 'LineWidth', 3);
        xlim([150, 650]);
        set(gca, 'YScale', 'log', 'color', 'none');
        title('(iv) Pre PTZ (WT)');

        % MUT enlarged plot.
        subplot(3, 2, 4);
        hold on;
        plot(x, y2, 'r-.', 'LineWidth', 1);
        plot(x(7:13), y2(7:13), 'r-', 'LineWidth', 3);
        xlim([150, 650]);
        set(gca, 'YScale', 'log', 'color', 'none');
        title('(v) Pre PTZ (scn1lab^{-/-})');
    end
end

% Add common labels for the entire figure.
han = axes(fig, 'visible', 'off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han, 'Correlation');
xlabel(han, 'Euclidean Distances (\mum)');

% 2. STATISTICAL COMPARISON OF FIT METRICS
% Create a new figure for the boxplots.
figure;

% Loop to plot both slope coefficient and RMSE.
for i = 1:2
    if i == 1
        data = slopecoef;
        ylabelnames = 'Power coefficient';
    else
        data = rmseVal;
        ylabelnames = 'RMSE';
    end

    % Select the data for the 'late PTZ' stage (kwid=3).
    subplot(1, 2, i);
    hold on;
    data = data(:, 3); % Change this value to compare other conditions.
    data1 = data(1:length(WT_fish));
    data2 = data(length(WT_fish) + 1:end);

    % Create boxplots with individual data points.
    boxplot(data1, 'Positions', 1, 'Symbol', '');
    boxplot(data2, 'Positions', 2, 'Symbol', '');
    scatter(1, data1, 'b.');
    scatter(2, data2, 'r.');
    yline(0, '--');

    % Plot the significance indicator.
    if h
        plot([1, 2], [max(data), max(data)], 'k');
        text(1.5, max(data) * 1.1, sig_star);
    end

    % Set axis limits and labels.
    if max(data) > 0
        ylim([min(data), max(data) * 1.2]);
    else
        ylim([min(data), 0]);
    end
    xlim([0.5, 2.5]);
    xticks([1, 2]);
    ylabel(ylabelnames);
    xticklabels({'WT', 'MUT'});
end

% Set common y-axis limits for the second plot.
ylim([0, 0.2]);

%% Panel C 

% This script compares the correlation-distance relationship between WT and
% SCN1 mutant fish across three different regional groupings:
% ipsilateral (side), contralateral (cross), and whole-brain.
% It normalizes the data, performs ANOVA, and plots the results.

%% 1. INITIALIZATION & DATA LOADING
clear all;
clc;

% Define the three data files to be analyzed.
files = {'side', 'cross', 'whole'};

% Initialize matrices to store statistical results.
varMat = [];
meanMat = [];
pMat = [];

% Loop through each file (regional grouping).
for fj = 1:length(files)
    % Load the data for the current file.
    load(['paper_' files{fj} '20.mat']);

    % Combine the last three distance bins (indices 7 to end) into a single bin.
    temp = [corrdisMat(1:6, :, :); nanmean(corrdisMat(7:end, :, :), 1)];
    corrdisMat = temp;

    % Normalize the data.
    % The normalization is performed for a specific condition (`kwid=9`),
    % which corresponds to the "late PTZ" stage in the previous script.
    kwid = 9;
    corrdisMatNorm = zeros(size(corrdisMat, 1) - 1, size(corrdisMat, 2));
    
    % Loop through each fish to apply the normalization.
    for fi = 1:size(corrdisMat, 2)
        % Normalization: Subtract the correlation at the smallest distance from all others.
        corrdisMatNorm(:, fi) = corrdisMat(2:end, fi, kwid) - nanmean(corrdisMat(1, fi, kwid));
    end

    % 2. STATISTICAL ANALYSIS (ANOVA)
    % Perform a two-way repeated measures ANOVA on the normalized data.
    statoutput = ranova2(corrdisMatNorm');

    % Loop through each distance bin to compute the difference in mean and variance.
    for i = 1:size(corrdisMatNorm, 1)
        % Extract the normalized data for the WT and mutant groups.
        data1 = reshape(corrdisMatNorm(i, 1:length(WT_fish), :), 1, []);
        data2 = reshape(corrdisMatNorm(i, length(WT_fish) + 1:end, :), 1, []);
        
        % Calculate and store the relative difference in variance and mean.
        varMat(i, fj) = (nanvar(data1) - nanvar(data2)) / (nanvar(data1) + nanvar(data2));
        meanMat(i, fj) = (nanmean(data1) - nanmean(data2)) / (nanmean(data1) + nanmean(data2));
    end
    
    % Store the p-values from the ANOVA for the current regional grouping.
    pMat(:, fj) = statoutput.pValue(2:2:end);
end

% 3. VISUALIZATION OF RESULTS
close all;
figure('position', [1000, 888.3, 453, 249.4]);
hold on;

% Create a colormap based on the `meanMat` values.
evensizemeaMat = [-meanMat, meanMat];
[~, colsidx] = sort(reshape(evensizemeaMat, 1, []));
[~, colsidx] = sort(colsidx);
cols = redblueRealTecplot(length(colsidx));
colMat = reshape(colsidx, size(evensizemeaMat));
colMat = colMat(:, 1:size(meanMat, 2));

% Plot the results.
for fj = 1:size(pMat, 2)
    for i = 1:size(pMat, 1)
        % Plot colored, filled circles representing the statistical results.
        % The size of the circle is proportional to the negative logarithm of the p-value.
        % The color represents the sign of the mean difference.
        scatter(i * 2 - 1, fj, -log(pMat(i, fj)) * 90, cols(colMat(i, fj), :), 'filled');
        
        % If the p-value is significant (< 0.05), draw a black circle around the point.
        if pMat(i, fj) < 0.05
            scatter(i * 2 - 1, fj, -log(pMat(i, fj)) * 90, [0, 0, 0], 'LineWidth', 2);
        end
    end
end

% 4. FINAL PLOT FORMATTING
% Set axis limits and ticks.
xlim([0, i * 2]);
ylim([0, fj + 1]);

% Set x-axis labels to represent the distance bins.
xticks(1:2:2 * size(corrdisMatNorm, 1));
xticklabels({'100', '200', '300', '400', '500', '>600'});
xlabel('Euclidean distance between ROIs (\mum)');

% Set y-axis labels to represent the regional groupings.
yticks(1:length(files));
yticklabels({sprintf('Ipsilateral hemisphere'), sprintf('Contralateral hemisphere'), sprintf('Whole brain \\newline(contra+ipsi)')});
ytickangle(30);

% Create a custom legend.
h = legend({'', '', '\mu_{WT}>\mu_{MUT}', '', '', '', '', '', '', '', 'significant', '\mu_{WT}<\mu_{MUT}', '', '', ''}, 'location', 'northoutside');

%% Panel E

% This script loads and processes pre-computed neural connectivity data
% to generate alluvial plots. The plots compare the mean absolute correlation
% differences between WT and SCN1 mutant fish for different brain regions.

% 1. INITIALIZATION AND DATA LOADING
clear all;
clc;
close all;

% Set a scaling factor for the plot flow thickness.
scaler = 42000;

% Define the stage of the experiment (1=Pre, 2=Early, 3=Late PTZ).
stg = 3;
if stg == 1
    stgrange = 1;
    load('chord_all_connections20_pre.mat');
elseif stg == 2
    stgrange = 2:8;
    load('chord_all_connections20_early.mat');
else
    stgrange = 9:20;
    load('chord_all_connections20_late.mat');
end

% Define the order of brain regions for plotting.
orderidx = [4, 1, 2, 3];

% 2. CONTRALATERAL CONNECTIONS (RIGHT SUBPLOT)
% Subplot for contralateral connections (between hemispheres).
subplot(1, 2, 2);

% Calculate the mean correlation difference between WT and MUT.
temp = nanmean(leftdataMatCross(:, :, :, stgrange), 4);
data = nanmean(temp(:, :, 1:length(WT_fish)), 3) - nanmean(temp(:, :, 1 + length(WT_fish):end), 3);

% Get color data from pre-computed tables.
% `outputctable` and `outputtable` likely contain statistical results
% (e.g., t-statistic or p-value) used for coloring.
lrid = 6;
colordata = (1 * (outputctable{lrid} > 0) - 1 * (outputctable{lrid} < 0)) .* outputtable{lrid};

% Rescale data and reorder based on `orderidx`.
data = abs(data(orderidx, orderidx)) * scaler * 2;
colordata = colordata(orderidx, orderidx);

% Define labels for the alluvial plot.
left_labels = {"Dien.(L)", "Mes.(L)", "Rho.(L)", "Tel.(L)"};
left_labels = left_labels(orderidx);
right_labels = {"Dien.(R)", "Mes.(R)", "Rho.(R)", "Tel.(R)"};
right_labels = right_labels(orderidx);

% Generate the alluvial plot for contralateral connections.
alluvialflowDiv(data, colordata, left_labels, right_labels, ['Contralateral'], false);
box off;

% 3. IPSILATERAL CONNECTIONS (LEFT SUBPLOT)
% Subplot for ipsilateral connections (within hemispheres).
subplot(1, 2, 1);

% Calculate the mean correlation difference for ipsilateral connections.
% Data from both left and right hemispheres are combined.
temp = nanmean(leftdataMatSame + rightdataMatSame, 4);
data = nanmean(temp(:, :, 1:length(WT_fish)), 3) - nanmean(temp(:, :, 1 + length(WT_fish):end), 3);

% Get color data from pre-computed tables.
lrid = 3;
colordata = (1 * (outputctable{lrid} > 0) - 1 * (outputctable{lrid} < 0)) .* outputtable{lrid};

% Rescale data and reorder.
data = abs(data(orderidx, orderidx)) * scaler;
colordata = colordata(orderidx, orderidx);

% Define labels.
left_labels = {"Dien.(L)", "Mes.(L)", "Rho.(L)", "Tel.(L)"};
left_labels = left_labels(orderidx);
right_labels = {"Dien.(L)", "Mes.(L)", "Rho.(L)", "Tel.(L)"};
right_labels = right_labels(orderidx);

% Generate the alluvial plot for ipsilateral connections.
alluvialflowDivSide(data, colordata, left_labels, right_labels, ['Ipsilateral'], false);
box off;
set(gca, 'xdir', 'reverse');
