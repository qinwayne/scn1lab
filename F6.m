%% Panel D and E
% 1. Initialization and Stage Configurations
clear all; clc; close all;
load('etagammarelationship_500.mat');

% Generate stage windows and labels ONCE for the entire script
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

num_WT = length(WT_fish);

% 500 simulations + 2D regression (Bubble Plot)
figure(3); hold on;

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
        p = pMat(kwid);
        if p < 0.05 && p >= 0.01
            txtMat{kwid} = '*';
        elseif p < 0.01 && p >= 0.001
            txtMat{kwid} = '**';
        elseif p < 0.001
            txtMat{kwid} = '***';
        end
    end

    % Color and size mapping
    evensizemeaMat = -[meanMat -meanMat];
    [~, colsidx] = sort(reshape(evensizemeaMat, 1, []));
    [~, colsidx] = sort(colsidx);

    cols = redblueRealTecplot(length(colsidx));
    colMat = reshape(colsidx, size(evensizemeaMat));
    colMat = colMat(:, 1:size(meanMat, 2));

    circlesize = -log(pMat)*40 + 0.0001;
    scatter(1:length(pMat), i, circlesize, cols(colMat,:), 'filled');
    
    circlesize(cellfun(@isempty, txtMat)) = nan;
    scatter(1:length(pMat), i, circlesize, [0 0 0], 'LineWidth', 2);
end

ylim([0 3]); yticks([1 2]); yticklabels({'|\gamma|','|\eta|'});
set(gcf, "Position", [986.3333, 567.0000, 560.0000, 181.3333]);
xticks(1:3:size(stages_ptz,1)+1);
xticklabels(labelList(1:3:end));

% Detailed change of parameters
figure('Renderer', 'painters');

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

    for kwid = 1:size(data, 2)
        data1 = wt_group(:, kwid);
        data2 = mut_group(:, kwid);
        
        h1 = scatter(kwid, data1, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.1, 'DisplayName', 'fish (WT)');
        h2 = scatter(kwid, data2, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.1, 'DisplayName', 'fish (MUT)');
        
        if kwid ~= 1
            set(h1, 'HandleVisibility', 'off');
            set(h2, 'HandleVisibility', 'off');
        end

        plotdata1(:, kwid) = mean(data1, 'omitnan');
        plotdata2(:, kwid) = mean(data2, 'omitnan');
    end

    statoutput = ranova2(data);
    
    skipwin = 1;
    count1 = 0; count2 = 0;
    yvals = ylim;
    
    for kwid = 1:skipwin:size(data, 2)-skipwin
        p = statoutput.pValue(kwid*2);
        if p < 0.05 
            p1 = patch([kwid kwid kwid+skipwin kwid+skipwin]-0.5, [yvals(1) yvals(2) yvals(2) yvals(1)], 'k', 'EdgeColor', 'none');
            
            if p > 0.01
                p1.FaceAlpha = 0.15;
                if count1 == 0, set(p1, 'DisplayName', 'p<0.05'); count1 = 1; else, set(p1, "HandleVisibility", "off"); end
            else
                p1.FaceAlpha = 0.4;
                if count2 == 0, set(p1, 'DisplayName', 'p<0.01'); count2 = 1; else, set(p1, "HandleVisibility", "off"); end
            end
        end
    end

    SEM1 = std(wt_group, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(wt_group)));
    SEM2 = std(mut_group, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(mut_group)));

    movwin = 1;
    x = movmean(1:size(data,2), movwin, "omitnan");
    y1 = movmean(mean(plotdata1, 1, 'omitnan'), movwin, "omitnan");
    y2 = movmean(mean(plotdata2, 1, 'omitnan'), movwin, "omitnan");
    
    plot(x, y1, 'b', 'DisplayName', 'WT median', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(x, y2, 'r', 'DisplayName', 'MUT median', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    x2 = [x, fliplr(x)];
    fill(x2, [y1 + SEM1, fliplr(y1 - SEM1)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'SEM (WT)', 'HandleVisibility', 'off');
    fill(x2, [y2 + SEM2, fliplr(y2 - SEM2)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'SEM (MUT)', 'HandleVisibility', 'off');

    ylabel(ylabname);
    if di == 2
        xticks(1:3:size(stages_ptz,1));
        xticklabels(labelList(1:3:end));
        xlabel('Stages (s)');
    else
        xticks([]);
    end
end
