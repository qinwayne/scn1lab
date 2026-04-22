%% Initial Setup
clear all; clc; close all;

% Group Definitions
WT_fish = {'07','14','17','19','20','24','25','28','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};
all_fish = [WT_fish, SCN1_fish];

% Parameters
p_th = 0.1;
winsz = 300;
slidsz = winsz/2;

%% Sliding Window Definitions
stages_ptz(1,:) = 100:100+winsz;
stages_ptz(2,:) = stages_ptz(1,:) + winsz;
i = 3;
while 700 + slidsz*(i-2) + winsz < 4350
    stages_ptz(i,:) = 700 + slidsz*(i-3) : 700 + slidsz*(i-3) + winsz;
    i = i + 1;
end
num_stages = size(stages_ptz, 1);

%% Load Brain Regions
load('..\Zbrain_Masks.mat'); % Contains Zbrain_Masks
RegionList = unique(Zbrain_Masks(:,1));
RegionList = RegionList([1 3 4 6]); % Subset of regions of interest

%% Connectivity Processing Loop
for fi = 1:length(all_fish)
    fprintf('Processing Fish %d/%d: %s\n', fi, length(all_fish), all_fish{fi});
    
    geno = 'WT'; if fi > length(WT_fish), geno = 'Hom'; end
    fish = all_fish{fi};
    load(['..\data\SCN1LabData\', geno, '\raw_fish_std_fmt_', fish, '.mat']);
    
    ROI_centroids = fish_ROI_centroids;
    midline = mean(ROI_centroids(:,2), 1);
    
    % Get Region indices for this specific fish
    PerBrainRegions = getPerBrainRegions(Zbrain_Masks, ROI_centroids, RegionList);
    regions = fieldnames(PerBrainRegions);
    num_regions = length(regions);

    for kwid = 1:num_stages
        % Calculate Correlation
        df_segment = fish_stim_trains{1,1}(:, stages_ptz(kwid,:));
        corrMat = corrcoef(df_segment');
        corrMat(isnan(corrMat)) = 0;
        corrMat = corrMat - diag(diag(corrMat));

        % Thresholding
        Wtgt = threshold_proportional(corrMat, p_th);
        [idx1, idx2] = find(Wtgt ~= 0);

        % Lateralization Checks
        cross_check = (ROI_centroids(idx1,2)-midline) .* (ROI_centroids(idx2,2)-midline) < 0;
        is_left     = ROI_centroids(idx1,2) < midline;
        is_right    = ROI_centroids(idx1,2) > midline;

        % Define filters for Whole, Left, and Right
        filters = { ...
            idx1(~cross_check), idx2(~cross_check), idx1(cross_check), idx2(cross_check); ... % Whole
            idx1(~cross_check & is_left), idx2(~cross_check & is_left), idx1(cross_check & is_left), idx2(cross_check & is_left); ... % Left
            idx1(~cross_check & is_right), idx2(~cross_check & is_right), idx1(cross_check & is_right), idx2(cross_check & is_right) ... % Right
        };

        for f_type = 1:3
            s1 = filters{f_type, 1}; s2 = filters{f_type, 2};
            c1 = filters{f_type, 3}; c2 = filters{f_type, 4};
            
            stemp = zeros(num_regions); ctemp = zeros(num_regions);
            
            for r1 = 1:num_regions
                r1_idx = PerBrainRegions.(regions{r1}).idx;
                
                % Within-side connections
                s_mask = ismember(s1, r1_idx);
                % Cross-side connections
                c_mask = ismember(c1, r1_idx);
                
                for r2 = 1:num_regions
                    r2_idx = PerBrainRegions.(regions{r2}).idx;
                    stemp(r1,r2) = sum(ismember(s2(s_mask), r2_idx)) / length(ROI_centroids)^2;
                    ctemp(r1,r2) = sum(ismember(c2(c_mask), r2_idx)) / length(ROI_centroids)^2;
                end
            end
            
            % Assign to global matrices
            if f_type == 1 % Whole
                dataMatSame(:,:,fi,kwid) = stemp; dataMatCross(:,:,fi,kwid) = ctemp;
            elseif f_type == 2 % Left
                leftdataMatSame(:,:,fi,kwid) = stemp; leftdataMatCross(:,:,fi,kwid) = ctemp;
            else % Right
                rightdataMatSame(:,:,fi,kwid) = stemp; rightdataMatCross(:,:,fi,kwid) = ctemp;
            end
        end
    end
end

%% Statistical Comparison
% Selecting stage window of interest (e.g., 9:20)
kwid_select = 9:min(20, num_stages);

% Perform ANOVA/T-test comparisons across different lateralization types
[outputtable{1}, outputctable{1}] = computeStats(leftdataMatCross(:,:,:,kwid_select));
[outputtable{2}, outputctable{2}] = computeStats(rightdataMatCross(:,:,:,kwid_select));
[outputtable{3}, outputctable{3}] = computeStats(leftdataMatSame(:,:,:,kwid_select) + rightdataMatSame(:,:,:,kwid_select));
[outputtable{4}, outputctable{4}] = computeStats(leftdataMatSame(:,:,:,kwid_select));
[outputtable{5}, outputctable{5}] = computeStats(rightdataMatSame(:,:,:,kwid_select));
[outputtable{6}, outputctable{6}] = computeStats(leftdataMatCross(:,:,:,kwid_select) + rightdataMatCross(:,:,:,kwid_select));

%% Save Output
save('chord_all_connections20_pre.mat', 'outputtable', 'outputctable', 'regions', ...
    'dataMatSame', 'dataMatCross', 'leftdataMatSame', 'leftdataMatCross', 'rightdataMatSame', 'rightdataMatCross');

%% Helper Functions
function [pTable, diffTable] = computeStats(data)
    num_r = size(data,1);
    pTable = zeros(num_r); diffTable = zeros(num_r);
    
    if size(data, 4) > 1
        for i = 1:num_r
            for j = 1:num_r
                input = squeeze(data(i,j,:,:));
                % Note: ranova2 is a custom function in your workspace
                [multicomp, aov] = ranova2(input, 4); 
                pTable(i,j) = 1 / aov.pValue(2); % Store as 1/p for chord plotting
                diffTable(i,j) = multicomp.Difference(2);
            end
        end
    else
        for i = 1:num_r
            for j = 1:num_r
                wt = data(i,j,1:17); mut = data(i,j,18:end);
                [~, p] = ttest2(wt, mut);
                pTable(i,j) = 1 / p;
                diffTable(i,j) = mean(mut) - mean(wt);
            end
        end
    end
end