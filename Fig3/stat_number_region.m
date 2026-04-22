%% Initial Setup
clc; clear; close all;

% Group Definitions
WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};
all_fish = [WT_fish, SCN1_fish];

base_path = 'D:\OneDrive - The University of Melbourne\Research\SCN1\data\SCN1LabData\region_roi_information\';

% Pre-allocate based on expected dimensions
% Note: Using variable names from original script logic
num_wt = length(WT_fish);
num_mut = length(SCN1_fish);
N = zeros(1, length(all_fish)); 

%% Process WT Group
fprintf('Processing WT Fish...\n');
for fi = 1:num_wt
    fish = WT_fish{fi};
    load([base_path, 'regioninfo_', fish, '.mat']); % Loads perRegion and fish_ROI_centroids
    
    N(fi) = size(fish_ROI_centroids, 1);
    
    % Extract ROI counts per subregion
    [region_counts, names] = getSubregionlength2(perRegion);
    
    % Initialize 3D matrices on first run
    if fi == 1
        WTidx = zeros(size(region_counts,1), size(region_counts,2), num_wt);
        WTnames = names;
    end
    WTidx(:,:,fi) = region_counts;
end

%% Process MUT (SCN1) Group
fprintf('Processing MUT Fish...\n');
for fi = 1:num_mut
    fish = SCN1_fish{fi};
    load([base_path, 'regioninfo_', fish, '.mat']);
    
    N(fi + num_wt) = size(fish_ROI_centroids, 1);
    
    [region_counts, names] = getSubregionlength2(perRegion);
    
    % Initialize 3D matrices on first run
    if fi == 1
        MUTidx = zeros(size(region_counts,1), size(region_counts,2), num_mut);
        MUTnames = names;
    end
    MUTidx(:,:,fi) = region_counts;
end

%% Save Essentials
% You can save as stat_number_Region.mat or stat_number_region_shortlist.mat
save('stat_number_Region.mat', 'WTidx', 'MUTidx', 'WTnames', 'MUTnames', 'N', 'WT_fish', 'SCN1_fish');
fprintf('Success: stat_number_Region.mat generated.\n');

%% Required Function
function [region_idx, region_name] = getSubregionlength2(perRegion)
    % Extracts the number of ROIs found in each predefined brain region
    regionList = {'Telencephalon','Diencephalon','Mesencephalon','Rhombencephalon'};
    
    % Determine max subregions to maintain consistent matrix size
    max_sub = 0;
    for r = 1:length(regionList)
        max_sub = max(max_sub, length(fieldnames(perRegion.(regionList{r}))));
    end
    
    region_idx = zeros(length(regionList), max_sub);
    region_name = cell(length(regionList), max_sub);
    
    for regi = 1:length(regionList)
        regionname = regionList{regi};
        subregionnames = fieldnames(perRegion.(regionname));
        for subi = 1:length(subregionnames)
            region_idx(regi,subi) = length(perRegion.(regionname).(subregionnames{subi}).idx);
            region_name{regi,subi} = [regionname ':' subregionnames{subi}];
        end
    end
end