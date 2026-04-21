clc; clear all
%% Panel B. Normalized Calcium Traces
% Grouping WT and SCN1 mutant fish identifiers
WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};
all_fish = [WT_fish, SCN1_fish];

% Prepare sliding windows
winSz = 300;
slidSz = winSz / 2;
stages_ptz(1, :) = 100:(100 + winSz);
allCond{1} = 'pre';

idx = 2;
while (1000 + slidSz * (idx - 1) + winSz) < 4200
    stages_ptz(idx, :) = (1000 + slidSz * (idx - 2)) : (1000 + slidSz * (idx - 2) + winSz);
    allCond{idx} = sprintf('stage%d', idx - 1);
    idx = idx + 1;
end

% Load specific fish index (matching the original loop setup for fi=20)
fi = 20; 
if fi > length(WT_fish)
    geno = 'Hom';
    subVal = 2;
else
    geno = 'WT';
    subVal = 1;
end

fishStr = all_fish{fi};
load('raw_fish_std_fmt_08.mat');

% 4. Raw Data Plot (PRODUCES PANEL B)
figure('Name', 'Panel B: Normalized Calcium Traces');
hold on;

window_range = 1900:4200;
% Create time vector and apply baseline offset
timeVec = (0:0.5:(size(fish_Suite2p_traces, 2) / 2 - 0.5)) - 300;
t_range = timeVec(window_range);

% Plot 7 staggered traces
for i = 1:7
    traceData = fish_Suite2p_traces(i * 2000, window_range);
    plot(t_range, normalize(traceData, 'range') + i, 'k');
end

xlim([min(t_range), max(t_range)]);
xticks(700:200:1700);
yticks([]);
xlabel('Time (s)');
ylabel('Normalised \DeltaF/F');

%% Panel C: Raster Plot, Active Neurons, and Tail Movement
% Clear previous plots
close all;

% Configuration
fishlist = [8, 20]; % e.g., fish ID 8 (WT) and 20 (Hom)
baseDir = '.\';
fsv = 199.94; % Frame sampling velocity
tailSeg = 3:5; % Tail segments for analysis

for fi = fishlist
    % Determine genotype, plotting color, and dynamic file paths
    if fi > length(WT_fish) % Assumes WT_fish is defined in Part 3
        geno = 'Hom'; 
        col = 'r'; % Red for SCN1 mutant
    else
        geno = 'WT';  
        col = 'b'; % Blue for WT
    end

    fishID = all_fish{fi};
    tailBasePath = fullfile(baseDir, sprintf('Fish_%s_base_beh.h5', fishID));
    tailPTZPath = fullfile(baseDir, sprintf('Fish_%s_5mM_beh.h5', fishID));

    % Initialize Figure
    figure('Position', [450, 230, 780, 420], 'Name', sprintf('Panel C - Fish %s (%s)', fishID, geno));

    %% 5A. Neural Data Loading & Binarization
    df_all = fish_stim_trains{1, 1};
    
    % Binarize (1 STD above mean) and calculate % active neurons
    threshold = mean(df_all, 'all') + std(df_all, 0, 'all');
    df_bin = df_all > threshold;
    actperc = (sum(df_bin, 1) / size(df_all, 1)) * 100; % Multiplied by 100 to match 0-100 figure axis

    %% 5B. Subplot 1: Raster Plot (Rows 1-5)
    subplot(7, 1, 1:5);
    imagesc(df_bin);
    colormap(flipud(gray(64))); % Black/White raster
    xline(600, '-.'); % PTZ addition marker
    xlim([1, 4200]);
    xticks([]);
    ylabel('Neurons');

    %% 5C. Subplot 2: % Active Neurons (Row 6)
    subplot(7, 1, 6);
    time_sec = 0.5:0.5:2100;
    bar(time_sec, actperc, 'FaceColor', col, 'EdgeColor', 'none');
    xline(300, '-.'); % 600 frames = 300 seconds
    xlim([1, 2100]);
    xticks([]);
    ylabel({'% Act.'; 'neuron'});

    %% 5D. Subplot 3: Tail Movement Trace (Row 7)
    subplot(7, 1, 7);
    hold on;
    
    if isfile(tailBasePath) && isfile(tailPTZPath)
        % Baseline phase
        tail_theta_base = h5read(tailBasePath, '/tail_theta');
        tail_theta_base = tail_theta_base - mean(tail_theta_base);
        tailBaseData = sum(abs(tail_theta_base(2000:end - 4001, tailSeg)), 2);
        tBase = (0:size(tailBaseData, 1) - 1) / fsv;
        plot(tBase, zscore(tailBaseData), 'k'); % Plotted in black to match figure

        % 5mM PTZ phase
        tail_theta_ptz = h5read(tailPTZPath, '/tail_theta');
        tailPTZData = sum(abs(tail_theta_ptz(1:end - 6000, tailSeg)), 2);
        tPTZ = (1:size(tailPTZData, 1)) / fsv + tBase(end);
        plot(tPTZ, zscore(tailPTZData), 'k'); 
        
        xlim([1, 2100]);
        ylim([-50, 50]);
    end
    ylabel({'Tail', 'angle^\circ'});

    %% 5E. X-Axis Formatting (Applied only to the bottom plot of the final fish)
    if fi == 20
        x_locs = 200:200:2100;
        x_labels = string(x_locs - 300); % Vectorized label creation
        x_labels(1) = "Pre PTZ";
        
        xticks(x_locs);
        xticklabels(x_labels);
        xlabel('Time (seconds post PTZ)');
    else
        xticks([]); % Hide X-axis for the upper plot (WT)
    end
end