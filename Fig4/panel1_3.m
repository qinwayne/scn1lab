%% Setup and Parameters
clc; clear; close all;

% Group Definitions
WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};
all_fish = [WT_fish, SCN1_fish];

% PTZ Stage Definitions
pre_ptz   = 100:400;
early_ptz = 600:2400;
late_ptz  = 2401:4200;
allcond   = {'pre', 'early', 'late'};
stages    = {pre_ptz, early_ptz, late_ptz};

% Analysis Parameters
distint = 30;           % Distance bin size
ROI_step = 4;           % Subsampling ROIs (1:4) to save memory/time
base_path = 'D:\OneDrive - The University of Melbourne\Research\SCN1\data\SCN1LabData\';
prop_path = 'D:\OneDrive - The University of Melbourne\Research\SCN1\paper\propagation\data\';

% Pre-allocate Result Matrices
num_fish = length(all_fish);
slopecoef = zeros(num_fish, 3);
rSqured   = zeros(num_fish, 3);
rmseVal   = zeros(num_fish, 3);
% corrdisMat will grow dynamically based on max distance found

%% Main Analysis Loop
for fi = 1:num_fish
    fprintf('Processing Fish %d/%d: %s\n', fi, num_fish, all_fish{fi});
    
    % 1. Load Data
    geno = ifthen(fi > length(WT_fish), 'Hom', 'WT');
    fish = all_fish{fi};
    
    % Load raw traces and propagation-derived centroids
    load(fullfile(base_path, geno, ['raw_fish_std_fmt_' fish '.mat']), 'fish_stim_trains');
    load(fullfile(prop_path, [fish '.mat']), 'raw_ROI_centroids');
    
    % 2. ROI Subsampling
    ROIbrainIdx = 1:ROI_step:size(raw_ROI_centroids, 1);
    centroids = raw_ROI_centroids(ROIbrainIdx, :);
    
    % 3. Distance Matrix Calculation
    distMat = pdist2(centroids, centroids);
    distMat(distMat == 0) = nan; 
    
    % Vectorize and sort distances once per fish
    [distVec, distInd] = sort(distMat(:));
    
    % 4. Condition Loop (Pre, Early, Late)
    for kwid = 1:3
        % Extract condition-specific traces
        df_all = fish_stim_trains{1, 1}(ROIbrainIdx, stages{kwid});
        
        % Correlation Calculation
        corrMat = corrcoef(df_all');
        corrMat(isnan(corrMat)) = 0;
        corrMat = corrMat - diag(diag(corrMat));
        
        % Vectorize and reorder correlations by distance
        corrVec = corrMat(:);
        ordCorrVec = corrVec(distInd);
        
        % 5. Data Cleaning
        nan_mask = isnan(distVec) | isnan(ordCorrVec);
        d_clean = distVec(~nan_mask);
        c_clean = ordCorrVec(~nan_mask);
        
        % 6. Distance Binning
        d_min = min(d_clean);
        d_max = max(d_clean);
        bins = d_min : distint : d_max - distint;
        
        for b = 1:length(bins)
            winidx = (d_clean > bins(b)) & (d_clean <= bins(b) + distint);
            if any(winidx)
                corrdisMat(b, fi, kwid) = mean(c_clean(winidx));
            else
                corrdisMat(b, fi, kwid) = nan;
            end
        end
        
        % 7. Power Law / Linear Fitting
        PLfits = getPLcoef(d_clean, c_clean);
        slopecoef(fi, kwid) = PLfits.Coefficients{2, 1};
        rSqured(fi, kwid)   = PLfits.Rsquared.Ordinary;
        rmseVal(fi, kwid)   = PLfits.RMSE;
    end
end

%% Save Final Variables
% Variables included: raw results, group lists, and stage definitions
save('panel1_3.mat', 'corrdisMat', 'slopecoef', 'rSqured', 'rmseVal', ...
     'WT_fish', 'SCN1_fish', 'all_fish', 'allcond', 'distint', 'ROI_step');

fprintf('Successfully generated panel1_3.mat\n');

%% Helper Functions
function val = ifthen(cond, a, b)
    if cond, val = a; else, val = b; end
end

function coef = getPLcoef(d, c)
    % Fits a linear model to correlation vs distance
    % Assuming d = distance, c = correlation
    coef = fitlm(d, c);
end