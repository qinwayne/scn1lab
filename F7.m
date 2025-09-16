%% Panel C and E
% This script analyzes and visualizes the results of a PLS-DA (Partial
% Least Squares Discriminant Analysis) model. It performs two main tasks:
% 1. Plots the decision boundary from a Support Vector Machine (SVM) on PLS-DA
%    scores.
% 2. Visualizes the contribution of different brain regions to the model's
%    class separation over time.

%% 1. Plot PLS-DA Scores with SVM Boundary
% This section iterates through timepoints to plot PLS-DA scores and a linear SVM boundary.
% The loop is intentionally kept to 1 for a static plot, as animating it is
% not standard for GitHub deployment unless the loop is for generating
% individual figures or a video.
for kwid = 1
    plot_svm_boundary(kwid);
    % drawnow;
    % pause(2);
end

%% 2. Analyze and Plot Region Contribution
% This section loads feature contribution data from PLS-DA and visualizes
% it as a heatmap.

select_var_path = 'geno'; 
subregion_path = '\subregiondata\'; 

% Specify the number of top features to analyze
top_features = 15000;

% Define region display indices, you can modify this to select specific regions
display_indices = 1:16; 

plot_region_contribution(select_var_path, subregion_path, top_features, display_indices);

function plot_svm_boundary(kwid)
% PLOT_SVM_BOUNDARY
% Reads PLS-DA score data, fits a linear SVM, and plots the data points
% with the SVM decision boundary.
%
% Input:
%   kwid: The identifier for the timepoint/file to load (e.g., 1, 2, ..., 20).

    % Load PLS-DA scores from CSV
    file_path = fullfile('geno', ['coef' num2str(kwid) '.csv']);
    if ~exist(file_path, 'file')
        error('File not found: %s', file_path);
    end
    data = readtable(file_path);
    
    % Prepare data for SVM
    X = [data.x, data.y];
    
    % Fit a linear Support Vector Machine model
    SVMModel = fitcsvm(X, data.group, 'KernelFunction', 'linear');
    
    % Get model parameters for plotting the decision boundary
    beta = SVMModel.Beta;
    b = SVMModel.Bias;
    
    % Create a new figure and plot
    figure('Renderer', 'painters', 'Position', [1.0000 1.0217 0.5600 0.3163]*1000);
    hold on;
    
    % Plot data points with different colors/styles for each group
    cols = [[0 0 1]; [1 0 0]];
    styles = {'o', '^'};
    gscatter(X(:,1), X(:,2), data.group, cols, styles, 8);
    
    % Plot the linear decision boundary
    x1_range = linspace(min(X(:,1)), max(X(:,1)), 100);
    x2_boundary = -(beta(1) / beta(2) * x1_range) - b / beta(2);
    plot(x1_range, x2_boundary, '-.');
    
    % Add labels and legend
    xlabel('PLS-DA Component 1');
    ylabel('PLS-DA Component 2');
    legend('WT', 'scn1lab^{-/-}', 'SVM Boundary', 'Location', 'best');
    
    hold off;
end

function plot_region_contribution(data_path, subregion_path, top_features, display_indices)
% PLOT_REGION_CONTRIBUTION
% Loads region contribution data from PLS-DA, processes it, and plots
% heatmaps to visualize region importance over time.
%
% Inputs:
%   data_path: The folder containing the 'selectVar' CSV files.
%   subregion_path: The path to the folder with subregion data for names.
%   top_features: The number of top features to consider for each timepoint.
%   display_indices: A vector of indices for regions to display on the y-axis.

    num_timepoints = 20;
    
    % Get region names from the subregion data folder
    all_folders = dir(subregion_path);
    region_folder = {all_folders(3:end).name};
    
    % Initialize matrix to store contribution counts
    sub_ids = zeros(num_timepoints, length(region_folder));
    
    % Loop through timepoints to calculate region contributions
    for kwid = 1:num_timepoints
        file_path = fullfile(data_path, ['selectVar' num2str(kwid) '.csv']);
        if ~exist(file_path, 'file')
            warning('File not found: %s. Skipping.', file_path);
            continue;
        end
        
        % Import the data with specified options
        opts = delimitedTextImportOptions("NumVariables", 3);
        opts.DataLines = [2, Inf];
        opts.Delimiter = ",";
        opts.VariableNames = ["name", "valuevar", "comp"];
        opts.VariableTypes = ["double", "double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        data = readtable(file_path, opts);
        
        % Calculate contribution of each region
        top_data_names = data.name(1:min(top_features, end));
        unique_regions = unique(ceil(top_data_names / 10000));
        
        for i = 1:length(unique_regions)
            region_id = unique_regions(i);
            sub_ids(kwid, region_id) = sum(ceil(top_data_names / 10000) == region_id);
        end
    end
    
    % Generate time labels for x-axis
    time_labels = create_time_labels(num_timepoints);
    
    %% Plot Heatmap of Contribution Over Time
    figure('Renderer', 'painters', 'Position', [1.0000 0.8177 0.7557 0.4200]*1000);
    imagesc(sub_ids(:, display_indices)' / top_features);
    colormap(flipud(hot));
    clim([0 0.3]);
    yticks(1:length(display_indices));
    yticklabels(region_folder(display_indices));
    xticks(1:3:num_timepoints);
    xticklabels(time_labels(1:3:end));
    xlabel('Time (s)');
    title('Region Contribution Over Time');
    
    %% Plot Heatmap of Overall Contribution
    figure('Renderer', 'painters', 'Position', [1.0000 0.5817 0.2483 0.4200]*1000);
    imagesc(sum(sub_ids(:, display_indices), 1)' / top_features / num_timepoints);
    colormap(flipud(hot));
    clim([0 0.3]);
    yticks(1:length(display_indices));
    yticklabels(region_folder(display_indices));
    title('Overall Region Contribution');
    
    % Plot the contribution for the first timepoint as a separate figure.
    % This is based on the original code's final figure.
    figure('Renderer', 'painters', 'Position', [1.0000 0.5817 0.2483 0.4200]*1000);
    imagesc(sub_ids(1, display_indices)' / top_features);
    colormap(flipud(hot));
    clim([0 0.3]);
    yticks(1:length(display_indices));
    yticklabels(region_folder(display_indices));
    title('Contribution at pre PTZ');
end

function label_list = create_time_labels(num_stages)
% Helper function to generate time labels for the x-axis.
    winsz = 300;
    slidsz = winsz/2;
    stages_ptz(1,:) = 100:100+winsz;
    i = 2;
    while 1000+slidsz*(i-1)+winsz < 4200 && i <= num_stages
        stages_ptz(i,:) = 1000+slidsz*(i-2):1000+slidsz*(i-2)+winsz;
        i = i+1;
    end
    
    label_list = cell(1, num_stages);
    for i = 1:size(stages_ptz, 1)
        val1 = stages_ptz(i,1)/2-300;
        label_list{i} = mat2str(val1);
    end
    label_list{1} = 'Pre PTZ';
end