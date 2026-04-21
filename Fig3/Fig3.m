%% Panel B
% This script compares seizure event density maps between WT and SCN1
% mutant fish. It performs a voxel-wise statistical analysis and visualizes
% the results as a series of 2D slices.

% 1. INITIALIZATION & DATA LOADING
clear all;
load("density_comparison2.mat");
load("stat_number_region.mat");

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
% colormap(flip(redblueRealTecplot));

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

% 1. INITIALIZE FIGURE
% Create a new figure with a specified size.
fig = figure('Position', [1000, 944.3, 560, 293.3]);

% 2. LOOP THROUGH BRAIN REGIONS AND PLOT DATA
for i = 1:4
    if i > 2
        subplot(4, 1, [1:2]);
    else
        subplot(4, 1, [3:4]);
        % Add a legend only to the bottom subplot to avoid redundancy.
        legend('Location', 'best');
    end
    
    hold on;
    xlim([0, 11]);
    
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
    
    scatter((i*2+1) - 1, WTratio(:, i), 'Marker', '.', 'MarkerEdgeColor', 'blue', 'HandleVisibility', 'off');
    scatter((i*2+1), MUTratio(:, i), 'Marker', '.', 'MarkerEdgeColor', 'red', 'HandleVisibility', 'off');
    
    % Create boxcharts to visualize the data distribution.
    if i == 1
        boxchart(i*2 + ones(length(WT_fish), 1) - 1, WTratio(:, i), 'BoxFaceColor', 'b', 'MarkerStyle', '.', 'MarkerColor', 'b', 'DisplayName', 'WT');
        boxchart(i*2 + ones(length(SCN1_fish), 1), MUTratio(:, i), 'BoxFaceColor', 'r', 'MarkerStyle', '.', 'MarkerColor', 'r', 'DisplayName', 'MUT');
    else
        boxchart(i*2 + ones(length(WT_fish), 1) - 1, WTratio(:, i), 'BoxFaceColor', 'b', 'MarkerStyle', '.', 'MarkerColor', 'b', 'HandleVisibility', 'off');
        boxchart(i*2 + ones(length(SCN1_fish), 1), MUTratio(:, i), 'BoxFaceColor', 'r', 'MarkerStyle', '.', 'MarkerColor', 'r', 'HandleVisibility', 'off');
    end
    
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
% Initialization
close all;
load("stat_number_region_shortlist.mat");

fig = figure;
hold on;
count = 0;
[collabels, RegionList] = getColours();
xlabelnames = {};

% Process Regions and Subregions
for regi = 1:size(WTidx, 1)
    for subi = 1:size(WTidx, 2)
        % Calculate normalized WT and MUT values
        nowt = squeeze(WTidx(regi, subi, :)) ./ N(1:length(WT_fish))';
        nomut = squeeze(MUTidx(regi, subi, :)) ./ N(1+length(WT_fish):end)';
        
        % Only proceed if there are no zero values in the arrays
        if ~any(nowt == 0) && ~any(nomut == 0)
            [~, p] = kstest2(nowt, nomut);

            % If statistically significant
            if p < 0.001 
                count = count + 1;
                
                % Calculate normalized mutant values relative to WT mean
                norm_mut = (nomut / mean(nowt, 'omitnan')) * 100;
                
                % Plot scatter and boxchart
                scatter(ones(length(SCN1_fish), 1) * count, norm_mut, ...
                    'Marker', '.', 'MarkerEdgeColor', collabels{regi}, 'HandleVisibility', 'off');
                
                h2 = boxchart(ones(length(SCN1_fish), 1) * count, norm_mut, ...
                    'BoxFaceColor', collabels{regi}, 'MarkerStyle', '.', ...
                    'MarkerColor', collabels{regi}, 'DisplayName', 'MUT', 'HandleVisibility', 'off');

                % Only show MUT in legend once
                if count == 1
                    set(h2, 'HandleVisibility', 'on');
                end
                
                yline(100, '-.', 'HandleVisibility', 'off');
                plot(2*count, mean(nomut, 'omitnan'), 'k', 'HandleVisibility', 'off');

                % Format and store labels
                tempname = strrep(MUTnames{regi, subi}, '_', ' ');
                tempname = strrep(tempname, 'cephalon', '.');
                xlabelnames{count} = tempname;
                
                % Store indices (used externally)
                subregidx(count, :) = [regi, subi]; 
            end
        end
    end

    % Process the 4 Main Regions
    count = count + 1;
    
    WTratio = zeros(1, length(WT_fish));
    MUTratio = zeros(1, length(all_fish) - length(WT_fish));
    
    % Calculate WT ratios
    for fi = 1:length(WT_fish)
        WTratio(fi) = max(WTidx(regi, :, fi), [], 2) / N(fi);
    end
    
    % Calculate MUT ratios
    for fi = 1:length(MUTratio)
        % Offset the 'N' index by the length of WT_fish
        n_idx = fi + length(WT_fish);
        MUTratio(fi) = max(MUTidx(regi, :, fi), [], 2) / N(n_idx);
    end
    
    norm_mut_ratio = (MUTratio / mean(WTratio, 'omitnan')) * 100;
    
    % Plot Main Region Data
    scatter(ones(length(SCN1_fish), 1) * count, norm_mut_ratio, ...
        'Marker', '.', 'MarkerEdgeColor', collabels{regi}, 'HandleVisibility', 'off');
        
    boxchart(ones(length(SCN1_fish), 1) * count, norm_mut_ratio, ...
        'BoxFaceColor', collabels{regi}, 'MarkerStyle', '.', ...
        'MarkerColor', collabels{regi}, 'HandleVisibility', 'off');
        
    xlabelnames{count} = RegionList{regi};
end

% Final Figure Formatting
xticks(1:(count + 0.5));
xlim([0.5, count + 0.5]);
xticklabels(xlabelnames);
ylabel('MUT ROI number change (%)');


%% Panel E
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
    filepath = fullfile('.\', ...
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

close all
clear all
load("catchup.mat")
data = meancorrfishlist;
close all
clc
figure('position',[1177         811         560         377],'Renderer','painters')
for kwid = 1:size(data,2)
    if ~any([5:11] == kwid ) 
        data1 = data(1:length(WT_fish),kwid);
        data2 = data(length(WT_fish)+1:end,kwid);%+4000; power increase to have the same PTZ baseline
    else
        data1 = nan;
        data2 = nan;
    end

    hold on
    h1 = scatter(kwid,data1,10,'b','filled','MarkerFaceAlpha',0.1,'DisplayName','fish (WT)');
    h2 = scatter(kwid,data2,10,'r','filled','MarkerFaceAlpha',0.1,'DisplayName','fish (MUT)');
    if kwid ~= 1
        set(h1,'HandleVisibility', 'off');
        set(h2,'HandleVisibility', 'off');
    end

    plotdata1(:,kwid) = mean(data1);
    plotdata2(:,kwid) = mean(data2);
end
legend([h1(1),h2(1)],'Location','best')
% use 2way anova instead of ttest
clc
[statoutput] = ranova2(data);

% Reshape your data if needed: 35 fish × 75 conditions
% Create a table with repeated measures
fishIDs = (1:size(meancorrfishlist,1))';
condNames = strcat("Cond", string(1:size(meancorrfishlist,2)));
T = array2table(meancorrfishlist, 'VariableNames', condNames);
T.Fish = fishIDs;

% Fit repeated measures model
rm = fitrm(T, sprintf('%s-%s ~ 1', condNames{1}, condNames{end}), 'WithinDesign', table(condNames', 'VariableNames', {'Condition'}));

% Run ANOVA
ranovaOutput = ranova(rm);

% Extract SS terms
SS_effect = ranovaOutput.SumSq(1);  % Adjust index as needed
SS_total = sum(ranovaOutput.SumSq);

eta2 = SS_effect / SS_total;
cohen_f = sqrt(eta2 / (1 - eta2));

fprintf('Cohen''s f: %.4f\n', cohen_f);
%
skipwin = 1;
count1 = 0;
count2 = 0;
for kwid = 1:skipwin:size(data,2)-skipwin

    data1 = reshape(data(1:length(WT_fish),kwid:kwid+skipwin-1),1,[]);
    data2 = reshape(data(length(WT_fish)+1:end,kwid:kwid+skipwin-1),1,[]);%+4000; power increase to have the same PTZ baseline

    p = statoutput.pValue(kwid*2);

    if p<0.05
        p1=patch([kwid kwid kwid+skipwin kwid+skipwin]-0.5, [0 0.9 0.9 0], 'k');

        if p <0.05 && p>0.01 
            p1.FaceAlpha = 0.15;
            if count1 == 0
                set(p1,'DisplayName','p<0.05');
                count1 = 1;
            else
                set(p1,"HandleVisibility","off")
            end
        elseif p<0.01 
            p1.FaceAlpha = 0.4;
            if count2 == 0
                set(p1,'DisplayName','p<0.01');
                count2 = 1;
            else
                set(p1,"HandleVisibility","off")
            end
        end
        p1.EdgeColor = 'none';
    end
end
% plot(1:4,data(1:length(WT_fish),1:4)','b','LineStyle',':')
% plot(12:75,data(1:length(WT_fish),12:end)','b','LineStyle',':')
% plot(1:4,data(1+length(WT_fish):end,1:4)','r','LineStyle',':')
% plot(12:75,data(1+length(WT_fish):end,12:end)','r','LineStyle',':')

SEM1 = std(data(1:length(WT_fish),:))/sqrt(length(WT_fish));    
SEM2 = std(data(1+length(WT_fish):end,:))/sqrt(length(SCN1_fish));   

% ylim([500 58000])
movwin = 3;
x = movmean(1:size(data,2),movwin,"omitnan");
y1 = movmean(mean(plotdata1,1),movwin,"omitnan");
y2 = movmean(mean(plotdata2,1),movwin,"omitnan");
plot(x(1:4),y1(1:4),'b','DisplayName','WT median','LineWidth',1,'HandleVisibility','off')
plot(x(1:4),y2(1:4),'r','DisplayName','MUT median','LineWidth',1,'HandleVisibility','off')
plot(x(12:end),y1(12:end),'b','DisplayName','WT median','LineWidth',1,'HandleVisibility','off')
plot(x(12:end),y2(12:end),'r','DisplayName','MUT median','LineWidth',1,'HandleVisibility','off')
x2 = [x, fliplr(x)];

curve1 = y1 + SEM1;
curve2 = y1 - SEM1;
inBetween = [curve1, fliplr(curve2)];
fill([x2(1:4) x2(end-3:end)], [inBetween(1:4) inBetween(end-3:end)], 'b','EdgeColor','none','FaceAlpha',0.3,'DisplayName','SEM (WT)','HandleVisibility','off');
fill(x2(12:end-11), inBetween(12:end-11), 'b','EdgeColor','none','FaceAlpha',0.3,'HandleVisibility','off');

curve1 = y2 + SEM2;
curve2 = y2 - SEM2;
inBetween = [curve1, fliplr(curve2)];
fill([x2(1:4) x2(end-3:end)], [inBetween(1:4) inBetween(end-3:end)], 'r','EdgeColor','none','FaceAlpha',0.3,'DisplayName','SEM (MUT)','HandleVisibility','off');
fill(x2(12:end-11), inBetween(12:end-11), 'r','EdgeColor','none','FaceAlpha',0.3,'HandleVisibility','off');

xlim([0 size(data,2)])
kwid = 1;
while (kwid)*slidewin+winsize+100 < 4200
    starttime = (kwid-1)*slidewin+100;
    endtime = (kwid)*slidewin+winsize+100;
    if starttime<600
        timelabels{kwid} = 'Pre PTZ';
    else
        timelabels{kwid} = [mat2str(starttime/2-300) '~' mat2str(endtime/2-300)];
    end
    kwid = kwid + 1;
end
xticks(3:10:size(data,2))
xticklabels(timelabels(3:10:size(data,2)));
ylabel('Mean Correlation')
xlabel('Time (s)')
legend('Location','best')


function [region_idx, region_name] = getSubregionlength(ROI_locations)
    load('D:\OneDrive - The University of Melbourne\Research\SCN1\region_shortlist_variables.mat');
    regionList = {'Telencephalon', 'Diencephalon', 'Mesencephalon', 'Rhombencephalon'};
    
    for regi = 1:length(regionList)
        regionname = regionList{regi};
        PerBrainRegions = getPerBrainRegionsAll(Zbrain_shortlist, ROI_locations, regionname);
        subregionnames = fieldnames(PerBrainRegions.(regionname));
        
        for subi = 1:length(subregionnames)
            region_idx(regi, subi) = length(PerBrainRegions.(regionname).(subregionnames{subi}).idx);
            region_name{regi, subi} = [regionname ':' subregionnames{subi}];
        end
    end
end

function [region_idx, region_name] = getSubregionlength2(perRegion)
    regionList = {'Telencephalon', 'Diencephalon', 'Mesencephalon', 'Rhombencephalon'};
    
    for regi = 1:length(regionList)
        regionname = regionList{regi};
        subregionnames = fieldnames(perRegion.(regionname));
        
        for subi = 1:length(subregionnames)
            region_idx(regi, subi) = length(perRegion.(regionname).(subregionnames{subi}).idx);
            region_name{regi, subi} = [regionname ':' subregionnames{subi}];
        end
    end
end

function vals = combineRegions(vals, regidx)
    for regi = 1:size(vals, 1)
        for subi = 1:size(regidx, 2)
            if ~isempty(regidx{regi, subi})
                vals(regi, regidx{regi, subi}(1), :) = sum(vals(regi, regidx{regi, subi}, :));
                vals(regi, regidx{regi, subi}(2:end), :) = 0;
            end
        end
    end
end

