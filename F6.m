% ANALYZE_NETWORK_PARAMETERS
%
% This script analyzes the 'eta' and 'gamma' network parameters from
% zebrafish brain data to compare wild-type (WT) and mutant (MUT) groups.
% It performs statistical tests (repeated-measures ANOVA)
% and visualizes the results as a significance plot and a detailed
% time-series plot.

%% 1. Setup and Data Loading
clear all;
clc;
close all;

% Load the main data file
load('etagammarelationship_500.mat', 'WT_fish', 'gammamat', 'etamat');

% Define the two parameters to analyze
parameters = {gammamat', -etamat'};
param_names = {'|\gamma|', '|\eta|'};

%% 2. Process and Plot Significance Circles
% This function performs the statistical tests and generates the plot.
plot_significance_circles(parameters, param_names, WT_fish);

%% 3. Plot Detailed Parameter Changes
% This function generates the detailed time-series subplot for a single parameter.
plot_parameter_details(-etamat', '|\eta|', WT_fish);


function plot_significance_circles(parameters, param_names, WT_fish)
% PLOT_SIGNIFICANCE_CIRCLES
% Analyzes and visualizes statistical significance using a circle plot.
%
% Inputs:
%   parameters: A cell array of matrices, where each matrix contains
%               parameter data for the two groups.
%   param_names: A cell array of strings for the plot's y-axis labels.
%   WT_fish: A cell array of WT fish IDs to determine group split.

    figure('Renderer', 'painters', 'Position', [986.3333 567.0000 560.0000 181.3333]);
    hold on;

    num_stages = size(parameters{1}, 2);
    num_groups = size(parameters{1}, 1);
    num_wt = length(WT_fish);
    load('p.mat')

    for i = 1:length(parameters)
        param_mat = parameters{i};
        param_mat(param_mat == 0) = NaN; % Replace zeros with NaN

        mean_diff = NaN(1, num_stages);

        for kwid = 1:num_stages
            data1 = param_mat(1:num_wt, kwid);
            data2 = param_mat(num_wt + 1:end, kwid);

            mean_diff(kwid) = nanmean(data1) - nanmean(data2);
        end

        % Prepare data for plotting
        even_size_mean_mat = -[mean_diff -mean_diff];
        [~, cols_idx] = sort(reshape(even_size_mean_mat, 1, []));
        [~, cols_idx] = sort(cols_idx);
        
        % Use a custom colormap for visualization
        cols = redblueRealTecplot(length(cols_idx));
        col_mat = reshape(cols_idx, size(even_size_mean_mat));
        col_mat = col_mat(:, 1:num_stages);

        % Calculate circle size from p-values
        circle_size = -log10(p_mat) * 40;
        circle_size(isinf(circle_size)) = 500; % Cap size for very small p-values
        
        % Plot filled circles
        scatter(1:num_stages, i, circle_size, cols(col_mat, :), 'filled');
        
        % Plot black borders for significant circles
        significant_idx = find(~cellfun(@isempty, sig_stars));
        scatter(significant_idx, repmat(i, size(significant_idx)), ...
            circle_size(significant_idx), [0 0 0], 'LineWidth', 2);
    end

    % Final plot adjustments
    ylim([0 3]);
    yticks([1 2]);
    yticklabels(param_names);
    xlabel('Time (s)');
    title('Statistical Significance of Network Parameters');
    
    % Add time labels to the x-axis
    label_list = calculate_time_labels(num_stages);
    xticks(1:3:num_stages);
    xticklabels(label_list(1:3:end));
end

function label_list = calculate_time_labels(num_stages)
% Helper function to generate time labels for the x-axis.
    winsz = 300;
    slidsz = winsz/2;
    stages_ptz(1, :) = 100 : 100 + winsz;
    i = 2;
    while 1000 + slidsz * (i - 1) + winsz < 4200 && i <= num_stages
        stages_ptz(i, :) = 1000 + slidsz * (i - 2) : 1000 + slidsz * (i - 2) + winsz;
        i = i + 1;
    end
    
    label_list = cell(1, num_stages);
    for i = 1:size(stages_ptz, 1)
        val1 = stages_ptz(i, 1) / 2 - 300;
        label_list{i} = mat2str(val1);
    end
    label_list{1} = 'Pre PTZ';
end

function plot_parameter_details(data_matrix, y_label_name, WT_fish)
% PLOT_PARAMETER_DETAILS
% Creates a subplot showing the mean and individual data points over time.
%
% Inputs:
%   data_matrix: A matrix of parameter data.
%   y_label_name: String for the y-axis label.
%   WT_fish: Cell array of WT fish IDs to determine group split.

    figure('Renderer', 'painters');
    
    num_wt = length(WT_fish);
    num_stages = size(data_matrix, 2);

    % Plot individual data points
    hold on;
    h1 = scatter(repmat(1:num_stages, num_wt, 1), data_matrix(1:num_wt, :), ...
        10, 'b', 'filled', 'MarkerFaceAlpha', 0.1, 'DisplayName', 'WT');
    h2 = scatter(repmat(1:num_stages, num_stages-num_wt, 1), data_matrix(num_wt+1:end, :), ...
        10, 'r', 'filled', 'MarkerFaceAlpha', 0.1, 'DisplayName', 'MUT');

    % Plot group means
    mean_wt = nanmean(data_matrix(1:num_wt, :), 1);
    mean_mut = nanmean(data_matrix(num_wt + 1:end, :), 1);
    
    plot(1:num_stages, mean_wt, 'b-', 'LineWidth', 2, 'DisplayName', 'Mean WT');
    plot(1:num_stages, mean_mut, 'r-', 'LineWidth', 2, 'DisplayName', 'Mean MUT');
    
    % Final plot adjustments
    xlabel('Time (s)');
    ylabel(y_label_name);
    title(['Detailed Change in ' y_label_name]);
    
    % Add legend
    legend([h1(1), h2(1)], 'WT', 'MUT');
    
    % Add time labels to the x-axis
    label_list = calculate_time_labels(num_stages);
    xticks(1:3:num_stages);
    xticklabels(label_list(1:3:end));
end