%% Initial Setup
clc; clear; close all;

% Group Definitions
WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};
all_fish = [WT_fish, SCN1_fish];

% Parameters
distint = 100;
distbins = 0:distint:600;
winsz = 300;
slidsz = winsz/2;

% Stage Generation
stages_ptz(1,:) = 100:100+winsz;
allcond{1} = 'pre';
i = 2;
while 1000 + slidsz*(i-1) + winsz < 4200
    stages_ptz(i,:) = 1000 + slidsz*(i-2) : 1000 + slidsz*(i-2) + winsz;
    allcond{i} = sprintf('stage%d', i-1);
    i = i + 1;
end

%% Main Processing Loop
% We run the analysis 3 times for the different "paper" outputs
analysis_types = {'whole', 'side', 'cross'};

for t = 1:length(analysis_types)
    type = analysis_types{t};
    fprintf('Starting analysis for: paper_%s20.mat\n', type);
    
    % Pre-allocate result matrix [dist_bins x fish x stages]
    corrdisMat = zeros(length(distbins), length(all_fish), length(allcond));
    
    for fi = 1:length(all_fish)
        fprintf('  Fish %d/%d: %s\n', fi, length(all_fish), all_fish{fi});
        
        % Load Data
        geno = ifthen(fi > length(WT_fish), 'Hom', 'WT');
        fish_path = ['D:\OneDrive - The University of Melbourne\Research\SCN1\data\SCN1LabData\', ...
                     geno, '\raw_fish_std_fmt_', all_fish{fi}, '.mat'];
        load(fish_path, 'fish_stim_trains', 'fish_ROI_centroids');
        
        ROI_centroids = fish_ROI_centroids;
        midline = mean(ROI_centroids(:,2));
        
        % Define Lateralization Masks
        is_left  = ROI_centroids(:,2) < midline;
        is_right = ROI_centroids(:,2) > midline;
        
        for kwid = 1:length(allcond)
            % Compute Correlations for the stage
            df_segment = fish_stim_trains{1,1}(:, stages_ptz(kwid,:));
            corr_full = corrcoef(df_segment');
            corr_full = corr_full - diag(diag(corr_full));
            
            % Compute Distances
            dist_full = pdist2(ROI_centroids, ROI_centroids);
            
            % Filter based on Analysis Type
            switch type
                case 'whole'
                    % No additional masking needed
                    C = corr_full; D = dist_full;
                case 'side'
                    % Only within-hemisphere connections
                    mask = (is_left * is_left') | (is_right * is_right');
                    C = corr_full; D = dist_full;
                    C(~mask) = nan; D(~mask) = nan;
                case 'cross'
                    % Only cross-hemisphere connections
                    C = corr_full(is_left, is_right);
                    D = dist_full(is_left, is_right);
            end
            
            % Vectorize and Clean
            C_vec = C(:); D_vec = D(:);
            valid = ~isnan(C_vec) & ~isnan(D_vec) & D_vec > 0;
            C_vec = C_vec(valid); D_vec = D_vec(valid);
            
            % Binning by Distance
            for di = 1:length(distbins)
                d_min = distbins(di);
                if d_min < 600
                    idx = D_vec > d_min & D_vec <= d_min + distint;
                else
                    idx = D_vec > d_min;
                end
                
                if any(idx)
                    corrdisMat(di, fi, kwid) = nanmedian(C_vec(idx));
                else
                    corrdisMat(di, fi, kwid) = nan;
                end
            end
        end
    end
    
    % Save specific file
    filename = sprintf('paper_%s20.mat', type);
    save(filename, 'corrdisMat', 'allcond', 'distbins', 'stages_ptz', 'WT_fish', 'SCN1_fish');
    fprintf('  Saved %s\n', filename);
end

%% Helper Functions
function val = ifthen(cond, a, b)
    if cond, val = a; else, val = b; end
end