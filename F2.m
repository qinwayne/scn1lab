% This script generates a comprehensive figure for a specific fish, combining
% a neural activity raster plot, the percentage of active neurons over time,
% and a corresponding tail movement trace. It is designed to compare a
% wild-type and a homozygous mutant fish.

%% 1. INITIALIZATION & DATA LOADING
% Clear previous plots.
close all;

% Define the fish IDs to be processed.
fishlist = [8, 20]; % Example: fish ID 8 (WT) and 20 (Hom).

% Loop through each specified fish.
for fi = fishlist
    % Create a new figure for each fish.
    fig = figure('Position', [447.6667, 232.3333, 782.0000, 420.0000]);
    
    % Determine the genotype and plot color based on the fish index.
    if fi > length(WT_fish)
        geno = 'Hom';
        col = 'red';
    else
        geno = 'WT';
        col = "blue";
    end

    % Load the neural activity data for the current fish.
    fish = all_fish{fi};
    filepath = fullfile('\The University of Melbourne\Research\SCN1\data\SCN1LabData', ...
                        geno, ['raw_fish_std_fmt_' fish '.mat']);
    load(filepath, 'fish_stim_trains');

    % Extract neural data and calculate binarized activity.
    df_all = fish_stim_trains{1, 1};
    % Binarize the data by setting values above a threshold (mean + 1 STD) to 1.
    df_bin = df_all > mean(df_all, 'all') + 1 * std(df_all, 0, 'all');
    
    % Calculate the percentage of active neurons at each time point.
    actperc = sum(df_bin, 1) / size(df_all, 1);

%% 2. PLOT NEURAL DATA
    % Subplot 1: Raster plot of binarized neural activity.
    subplot(7, 1, 1:5);
    imagesc(df_bin); % Display the binarized activity matrix.
    colormap(flipud(gray(64))); % Use a grayscale colormap.
    xline(600, '-.'); % Draw a vertical dashed line at the 600th frame.
    xlim([1, 4200]);
    xticks([]); % Hide x-axis ticks.
    ylabel('# Neuron');

    % Subplot 2: Bar plot of the percentage of active neurons over time.
    subplot(7, 1, 6);
    % The x-axis is scaled to represent time in seconds.
    bar(0.5:0.5:2100, actperc, 'FaceColor', col, 'EdgeColor', 'none');
    xlim([1, 4200] / 2);
    ylabel({'% Active'; 'neuron'});
    xticks([]); % Hide x-axis ticks.
    xline(600 / 2, '-.'); % Draw a vertical dashed line at the 300-second mark.

%% 3. PLOT TAIL MOVEMENT DATA
    % Subplot 3: Tail movement trace.
    subplot(7, 1, 7);
    hold on;
    
    % Define file paths and parameters for tail movement data.
    fsv = 199.94; % Frame sampling velocity.
    tailseg = 3:5; % Use points 3 to 5 of the tail segments for analysis.
    filename1 = fullfile('\The University of Melbourne\Research\SCN1\data\SCN1LabData\tail', ...
                         ['Fish_' all_fish{fi} '_base_beh.h5']);
    filename2 = fullfile('\The University of Melbourne\Research\SCN1\data\SCN1LabData\tail', ...
                         ['Fish_' all_fish{fi} '_5mM_beh.h5']);
    
    % Check if the data files exist before proceeding.
    if isfile(filename1) && isfile(filename2)
        extraframe = 6000;
        
        % Load and process baseline tail movement data.
        tail_theta = h5read(filename1, '/tail_theta');
        tail_theta = tail_theta - mean(tail_theta);
        data1 = sum(abs(tail_theta(2000:end - 4001, tailseg)), 2);
        tbase = (0:size(data1, 1) - 1) / fsv; % Time vector in seconds.
        plot(tbase, zscore(data1), 'g', 'HandleVisibility', 'off');

        % Load and process 5mM PTZ tail movement data.
        tail_theta = h5read(filename2, '/tail_theta');
        data2 = sum(abs(tail_theta(1:end - extraframe, tailseg)), 2);
        t = (1:size(data2, 1)) / fsv + tbase(end); % Time vector for PTZ phase.
        
        % Plot the tail movement trace.
        plot(t, zscore(data2), 'g', 'DisplayName', 'tail movement');
        xlim([1, 4200] / 2);
        ylim([-50, 50]);
    end

    % Format X-axis labels for the second fish (fi=20) for clarity.
    if fi == 20
        xlabellocs = 200:200:2100;
        xlabelnames = cell(1, length(xlabellocs));
        for i = 1:length(xlabellocs)
            xlabelnames{i} = mat2str(xlabellocs(i) - 300);
        end
        xlabelnames{1} = 'Pre PTZ';
        xticks(xlabellocs);
        xticklabels(xlabelnames);
        ylabel({'Tail', 'angle °'});
    end

%% 4. FINAL FIGURE FORMATTING
    % Add a common X-axis label for the entire figure.
    han = axes(fig, 'visible', 'off');
    han.Title.Visible = 'on';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    xlabel(han, 'Time (second)');
end