%% Initial Setup
clc; clear; close all;

% Group Definitions
WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};
all_fish = [WT_fish, SCN1_fish];
input_fish = all_fish;

%% Window and Stage Definition
winsz = 300;
slidsz = winsz/2;

% Initialize pre-PTZ stage
stages_ptz(1,:) = 100:100+winsz;
allcond{1} = 'pre';

% Generate PTZ stages
i = 2;
while 1000 + slidsz*(i-1) + winsz < 4200
    stages_ptz(i,:) = 1000 + slidsz*(i-2) : 1000 + slidsz*(i-2) + winsz;
    allcond{i} = sprintf('stage%d', i-1);
    i = i + 1;
end

%% Connectivity and Power Analysis
pwrfishlist = [];
meancorrfishlist = [];
slidewin = 50;
winsize = 300; 

for fi = 1:length(all_fish)
    fprintf('Processing Fish %d/%d: %s\n', fi, length(all_fish), all_fish{fi});
    
    % Determine Genotype path
    geno = 'WT';
    if fi > length(WT_fish), geno = 'Hom'; end
    
    fish = input_fish{fi};
    load(['D:\Wei\SCN1LabData\' geno '\raw_fish_std_fmt_' fish '.mat']);

    kwid = 1;
    % Loop through time segments
    while (kwid)*slidewin + winsize + 100 < size(fish_stim_trains{1,1}, 2)
        
        % Extract segment
        start_idx = (kwid-1)*slidewin + 100;
        end_idx   = (kwid)*slidewin + winsize + 100;
        df_all    = fish_stim_trains{1,1}(:, start_idx:end_idx);
        
        % Calculate Correlation Matrix
        corrMat = corrcoef(df_all'); 
        corrMat(isnan(corrMat)) = 0;
        corrMat = corrMat - diag(diag(corrMat)); % Remove self-correlation

        % Calculate Signal Power (Sum of squares)
        pwrfishlist(fi, kwid) = sum(abs(df_all).^2, 'all');

        % Calculate Mean Correlation
        meancorrfishlist(fi, kwid) = nanmean(corrMat, 'all');

        kwid = kwid + 1;
    end
end

%% Save All Variables
save('catchup.mat', '-v7.3');
fprintf('Success: catchup.mat generated with all variables.\n');