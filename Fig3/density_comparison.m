%% Initial Setup
clc; clear; close all;

% Group Definitions
WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};
all_fish = [WT_fish, SCN1_fish];

% Parameters
cubesize = 6;
num_planes = 24; % Adjusted based on your loop 0:24
base_path = '..\data\SCN1LabData\suite2p\';

% Pre-allocate cell arrays
fish_ROI_centroids = cell(length(all_fish), 1);
reszden = cell(length(all_fish), 1);
N = cell(length(all_fish), 1);

%% Part 1: Extract ROI Centroids from Suite2p
fprintf('Extracting ROI centroids...\n');

for fi = 1:length(all_fish)
    fprintf('Processing Fish %d/%d: %s\n', fi, length(all_fish), all_fish{fi});
    ROI_centroids = [];
    
    % Find the specific suite2p folder for this fish
    fishfolder = dir([base_path, 'suite2p_*fish', all_fish{fi}, '_*']);
    
    if isempty(fishfolder)
        warning('Folder not found for fish %s', all_fish{fi});
        continue;
    end
    
    % Loop through each imaging plane
    for i = 0:num_planes
        plane_file = fullfile(base_path, fishfolder.name, ['plane', num2str(i)], 'Fall.mat');
        
        if isfile(plane_file)
            load(plane_file, 'stat', 'iscell');
            cellidx = find(iscell(:, 1));
            
            % Calculate centroids for valid cells
            for ci = 1:length(cellidx)
                curr_stat = stat{cellidx(ci)};
                ROI_centroids = [ROI_centroids; [mean(curr_stat.xpix), mean(curr_stat.ypix), i]];
            end
        end
    end
    fish_ROI_centroids{fi} = ROI_centroids;
end

%% Part 2: Generate Density Distributions
fprintf('Generating density distributions...\n');

targetSize = []; % To be defined by the first fish processed

for fi = 1:length(all_fish)
    ROI_centroids = fish_ROI_centroids{fi};
    if isempty(ROI_centroids), continue; end
    
    % Use original indices or random permutation as per your logic
    randMat = randperm(size(ROI_centroids, 1));
    randLoc = ROI_centroids(randMat, :);
    
    % Calculate 3D density
    [N{fi}, pix_den_dis] = getROIsDensity(randLoc, cubesize);
    
    % Set target size for resampling based on the first valid fish
    if isempty(targetSize)
        targetSize = size(pix_den_dis);
    end
    
    % Resize 3D volume to match the target spatial dimensions
    reszden{fi} = imresize3(pix_den_dis, targetSize);
end

%% Save Essentials
save('density_comparison.mat', 'reszden', 'fish_ROI_centroids', 'targetSize', 'WT_fish', 'SCN1_fish', 'all_fish', '-v7.3');
fprintf('Success: density_comparison.mat generated.\n');

%% Required Function
function [N, avgMat] = getROIsDensity(data, cubesize)
    % N: total number of ROIs
    % avgMat: 3D matrix containing ROI counts per cube
    N = size(data, 1);
    
    xRange = min(data(:,1)):cubesize:max(data(:,1));
    yRange = min(data(:,2)):cubesize:max(data(:,2));
    zRange = min(data(:,3)):cubesize:max(data(:,3));
    
    avgMat = zeros(length(xRange), length(yRange), length(zRange));
    
    for xi = 1:length(xRange)-1
        xidx = data(:,1) >= xRange(xi) & data(:,1) < xRange(xi+1);
        for yi = 1:length(yRange)-1
            yidx = data(:,2) >= yRange(yi) & data(:,2) < yRange(yi+1);
            for zi = 1:length(zRange)-1
                zidx = data(:,3) >= zRange(zi) & data(:,3) < zRange(zi+1);
                
                roiCount = sum(xidx & yidx & zidx);
                if roiCount > 0
                    avgMat(xi, yi, zi) = roiCount;
                end
            end
        end
    end
end