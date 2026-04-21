%% 1. Initialization and Stage Configurations
clear all; clc; close all;

WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};
all_fish = [WT_fish, SCN1_fish];
num_WT = length(WT_fish);

% Generate stage windows and X-axis labels ONCE
winsz = 300; 
slidsz = winsz/2;
stages_ptz(1,:) = 100:100+winsz;
idx = 2;
while 1000 + slidsz*(idx-1) + winsz < 4200
    stages_ptz(idx,:) = 1000 + slidsz*(idx-2) : 1000 + slidsz*(idx-2) + winsz;
    idx = idx + 1;
end

labelList = cell(1, size(stages_ptz, 1));
for i = 1:size(stages_ptz, 1)
    labelList{i} = mat2str(stages_ptz(i,1)/2 - 300);
end
labelList{1} = 'Pre PTZ';

%% 2. Load Data for a Specific Region
allfolders = dir('subregiondata');
fileid = 3; % Change this index to loop through Pallium, Habenula, etc.
regionfolder = fullfile('subregiondata', allfolders(fileid).name);
addpath(fullfile(pwd, regionfolder));

etamat = zeros(20, length(all_fish));
gammamat = zeros(20, length(all_fish));

for fi = 1:length(all_fish)
    EbAll = [];
    for i = 1:30
        try
            fish = all_fish{fi};
            load(fullfile(regionfolder, [mat2str(i) '_' fish 'Eb20.mat']), 'Eb');
            load(fullfile(regionfolder, [mat2str(i) '_' fish '.mat']), 'params');
            EbAll(:,:,i) = squeeze(Eb(:,:,2)); 
        catch
            continue; % Skip missing files gracefully
        end
    end
    
    if ~isempty(EbAll)
        EbAll(EbAll==0) = nan;
        EbAvg = mean(EbAll, 3, 'omitnan');
        
        if ~isnan(EbAvg(1,1))
            for kwid = 1:size(EbAvg,2)
                % Find the indices of the 5 minimum unique values
                [unq_vals, ~] = unique(EbAvg(:,kwid));
                mink_vals = unq_vals(1:min(5, length(unq_vals)));
                Ebidx = find(ismember(EbAvg(:,kwid), mink_vals));
                
                etamat(kwid,fi) = mean(params(Ebidx,1));
                gammamat(kwid,fi) = mean(params(Ebidx,2));
            end
        end
    end
end

%% 3. PANEL A: Line Trace Plot with SEM and Patches (e.g., Pallium)
figure('Position', [100, 100, 600, 400], 'Renderer', 'painters');
sgtitle(['Region: ' allfolders(fileid).name], 'Interpreter', 'none');

for di = 1:2
    if di == 1
        data = -etamat'; ylabname = '|\eta|';
    else
        data = gammamat'; ylabname = '|\gamma|';
    end
    
    subplot(2, 1, di); hold on;
    
    wt_group = data(1:num_WT, :);
    mut_group = data(num_WT+1:end, :);
    
    plotdata1 = zeros(size(wt_group));
    plotdata2 = zeros(size(mut_group));

    % Scatter raw data
    for kwid = 1:size(data, 2)
        data1 = wt_group(:, kwid);
        data2 = mut_group(:, kwid);
        
        h1 = scatter(kwid, data1, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.1, 'DisplayName', 'WT');
        h2 = scatter(kwid, data2, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.1, 'DisplayName', 'scn1lab^{-/-}');
        
        if kwid ~= 1
            set(h1, 'HandleVisibility', 'off');
            set(h2, 'HandleVisibility', 'off');
        end

        plotdata1(:, kwid) = mean(data1, 'omitnan');
        plotdata2(:, kwid) = mean(data2, 'omitnan');
    end

    % Statistics & Patches
    statoutput = ranova2(data);
    yvals = ylim;
    skipwin = 1;
    
    for kwid = 1:skipwin:size(data, 2)-skipwin
        p = statoutput.pValue(kwid*2);
        if p < 0.05 
            p1 = patch([kwid kwid kwid+skipwin kwid+skipwin]-0.5, [yvals(1) yvals(2) yvals(2) yvals(1)], 'k', 'EdgeColor', 'none');
            if p > 0.01
                p1.FaceAlpha = 0.15; % Lighter gray for p<0.05
            else
                p1.FaceAlpha = 0.40; % Darker gray for p<0.01
            end
            set(p1, 'HandleVisibility', 'off');
        end
    end

    % SEM and Trendlines
    SEM1 = std(wt_group, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(wt_group)));
    SEM2 = std(mut_group, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(mut_group)));

    movwin = 1;
    x = movmean(1:size(data,2), movwin, "omitnan");
    y1 = movmean(mean(plotdata1, 1, 'omitnan'), movwin, "omitnan");
    y2 = movmean(mean(plotdata2, 1, 'omitnan'), movwin, "omitnan");
    
    plot(x, y1, 'b', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(x, y2, 'r', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    x2 = [x, fliplr(x)];
    fill(x2, [y1 + SEM1, fliplr(y1 - SEM1)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
    fill(x2, [y2 + SEM2, fliplr(y2 - SEM2)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');

    ylabel(ylabname);
    xlim([0.5, size(data,2)+0.5]);
    
    if di == 2
        xticks(1:3:size(stages_ptz,1));
        xticklabels(labelList(1:3:end));
        xlabel('Time (seconds post PTZ)');
    else
        xticks([]);
    end
end

% 4. PANEL B: Bubble Plot (Best Parameters Comparison)
figure('Position', [750, 100, 600, 200], 'Renderer', 'painters');
hold on;
title(['Region: ' allfolders(fileid).name], 'Interpreter', 'none');

for i = 1:2
    if i == 1
        paramat = gammamat';
    else
        paramat = -etamat';
    end
    paramat(paramat == 0) = nan;
    
    wt_data = paramat(1:num_WT, :);
    mut_data = paramat(num_WT+1:end, :);

    % Statistics
    [~, p1] = ttest2(wt_data(:,1), mut_data(:,1));
    pMat2 = ranova2(paramat(:,2:8)).pValue(2:2:end);
    pMat3 = ranova2(paramat(:,9:end)).pValue(2:2:end);
    pMat = [p1; pMat2; pMat3];

    meanMat = mean(wt_data, 1, 'omitnan') - mean(mut_data, 1, 'omitnan');
    txtMat = cell(1, size(paramat, 2));

    for kwid = 1:size(paramat, 2)
        if pMat(kwid) < 0.05
            txtMat{kwid} = '*';
        end
    end

    % Color mapping
    evensizemeaMat = -[meanMat -meanMat];
    [~, colsidx] = sort(reshape(evensizemeaMat, 1, []));
    [~, colsidx] = sort(colsidx);

    cols = redblueRealTecplot(length(colsidx)); 
    colMat = reshape(colsidx, size(evensizemeaMat));
    colMat = colMat(:, 1:size(meanMat, 2));

    % Plot Bubbles
    circlesize = -log(pMat)*40 + 0.0001;
    scatter(1:length(pMat), i, circlesize, cols(colMat,:), 'filled');
    
    % Outline significant bubbles
    circlesize(cellfun(@isempty, txtMat)) = nan;
    scatter(1:length(pMat), i, circlesize, [0 0 0], 'LineWidth', 1.5);
end

ylim([0.5 2.5]); 
yticks([1 2]); 
yticklabels({'|\gamma|','|\eta|'});
xticks(1:3:size(stages_ptz,1)+1);
xticklabels(labelList(1:3:end));
xlabel('Time (seconds post PTZ)');
xlim([0.5, length(pMat)+0.5]);

%% Panel C  D  E and F
% 1. Initialization and Time/Stage Labels
clear variables; close all; clc;

% Define stage windows and X-axis labels (Calculated once)
winsz = 300;
slidsz = winsz / 2;
stages_ptz(1, :) = 100:(100 + winsz);

idx = 2;
while 1000 + slidsz * (idx - 1) + winsz < 4200
    stages_ptz(idx, :) = 1000 + slidsz * (idx - 2) : 1000 + slidsz * (idx - 2) + winsz;
    idx = idx + 1;
end

labelList = cell(1, size(stages_ptz, 1));
for i = 1:size(stages_ptz, 1)
    labelList{i} = mat2str(stages_ptz(i, 1) / 2 - 300);
end
labelList{1} = 'Pre PTZ';

x_tick_locs = 1:3:size(stages_ptz, 1) + 1;
x_tick_labels = labelList(1:3:end);

% 2. PANELS C & D: PLS-DA Scatter Plots with SVM Boundary
% Loop through components (or specify a specific kwid for the plot)
kwid = 1; % Adjust to the specific stage you want to plot

% For Panel C (Genotype): Use 'geno\coefX.csv'
% For Panel D (Seizures): Change to your seizure data path (e.g., 'seizure\coefX.csv')
data = readtable(['D:\OneDrive - The University of Melbourne\Research\SCN1\PLSDA\geno\coef' mat2str(kwid) '.csv']); 

X = [data.x, data.y];
SVMModel = fitcsvm(X, data.group);

% Extract Linear predictor coefficients and Bias
beta = SVMModel.Beta; 
b = SVMModel.Bias; 

figure('Position', [100, 100, 560, 316], 'Renderer', 'painters');
hold on;

% Scatter plot for groups
cols = [0 0 1; 1 0 0]; % Blue (WT), Red (scn1lab-/-)
% Note: For Panel D, change cols to Green/Purple and styles to '*','.'
styles = {'o', '^'}; 
gscatter(X(:,1), X(:,2), data.group, cols, styles, 8);

% Calculate and plot the SVM Boundary Line
X1 = linspace(min(X(:,1)), max(X(:,1)), 100);
X2 = -(beta(1) / beta(2) * X1) - b / beta(2);
plot(X1, X2, '-.k', 'DisplayName', 'SVM');

legend({'WT', 'scn1lab^{-/-}', 'SVM'}, 'Location', 'best');
xlabel('PLS-DA Comp 1');
ylabel('PLS-DA Comp 2');
title(sprintf('PLS-DA Component Analysis (Stage %d)', kwid));
hold off;

% 3. PANEL E: Regional Classification Contribution Heatmap
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["name", "valuevar", "comp"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "name", "TrimNonNumeric", true);
opts = setvaropts(opts, "name", "ThousandsSeparator", ",");

subids = zeros(20, 16);
topfeatures = 15000;

for k = 1:20
    % Update path to match your local directory
    filepath = ['D:\OneDrive - The University of Melbourne\Research\SCN1\PLSDA\geno\selectVar' mat2str(k) '.csv'];
    if isfile(filepath)
        data_var = readtable(filepath, opts);
        temp = unique(ceil(data_var.name(1:topfeatures) / 10000));
        for i = 1:length(temp)
            subids(k, temp(i)) = sum(ceil(data_var.name(1:topfeatures) / 10000) == temp(i));
        end
    end
end

% Get Region Names
allfolders = dir('..\paper\Akarca_2023\subregiondata\');
regionfolder = {allfolders(3:end).name};
displayidx = 1:16; % Specific regions to display

figure('Position', [100, 500, 755, 420], 'Name', 'Panel E');
imagesc(subids(:, displayidx)' / topfeatures);
colormap(flipud(hot));
clim([0 0.33]);

% Formatting
yticks(1:16);
yticklabels(regionfolder);
xticks(x_tick_locs);
xticklabels(x_tick_labels);
xlabel('Time (seconds post PTZ)');

% Colorbar
cbar = colorbar;
cbar.Label.String = 'classification contribution';
cbar.Ticks = [0, 0.30];
cbar.TickLabels = {'0%', '30%'};

% 4. PANEL F: Overall Summaries (i, ii, iii, iv)
figure('Position', [900, 500, 250, 420], 'Name', 'Panel F');
colormap(flipud(hot));

% i: Genotype (pre PTZ)
subplot(1, 4, 1);
imagesc(sum(subids(1, displayidx), 1)' / topfeatures);
clim([0 0.33]);
yticks(1:16);
yticklabels(regionfolder);
xticks([]);
title('i');

% ii: Number of seizures (pre PTZ)
% (Replace `subids` with your seizure model output matrix if available)
subplot(1, 4, 2);
imagesc(sum(subids(1, displayidx), 1)' / topfeatures); % Placeholder using geno data
clim([0 0.33]);
yticks([]); xticks([]);
title('ii');

% iii: Genotype (overall)
subplot(1, 4, 3);
imagesc(sum(subids(:, displayidx), 1)' / topfeatures / 20);
clim([0 0.33]);
yticks([]); xticks([]);
title('iii');

% iv: Number of seizures (overall)
% (Replace `subids` with your seizure model output matrix if available)
subplot(1, 4, 4);
imagesc(sum(subids(:, displayidx), 1)' / topfeatures / 20); % Placeholder using geno data
clim([0 0.33]);
yticks([]); xticks([]);
title('iv');