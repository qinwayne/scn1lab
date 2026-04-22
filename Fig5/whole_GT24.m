%% Initial Setup & Parameters
clc; clear; close all;

% Group Definitions
WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};
all_fish = [WT_fish, SCN1_fish];

% Constants
swsize = 2000;      % Subsampling size for complex GT metrics
prop_threshold = 0.05; 
winsz = 300;
slidsz = winsz/2;

%% Window Generation (Stages)
stages_ptz(1,:) = 100:100+winsz;
stages_ptz(2,:) = stages_ptz(1,:) + winsz;
allcond = {'pre1', 'pre2'};

i = 3;
while 700 + slidsz*(i-2) + winsz < 4350
    stages_ptz(i,:) = 700 + slidsz*(i-3) : 700 + slidsz*(i-3) + winsz;
    allcond{i} = sprintf('stage%d', i-2);
    i = i + 1;
end

% Pre-allocate structure for data integrity
num_fish = length(all_fish);
num_cond = length(allcond);
gt = struct();

%% Main Processing Loop
for fi = 1:num_fish
    fprintf('\nProcessing Fish %d/%d: %s\n', fi, num_fish, all_fish{fi});
    
    % Determine Genotype and Path
    geno = ifthen(fi > length(WT_fish), 'Hom', 'WT');
    load(['..\..\SCN1\data\SCN1LabData\', geno, '\raw_fish_std_fmt_', all_fish{fi}, '.mat']);
    
    ROI_centroids = fish_ROI_centroids;
    N(fi) = size(ROI_centroids, 1);
    
    for kwid = 1:num_cond
        fprintf('  Calculating Stage: %s... ', allcond{kwid});
        
        % Correlation Analysis
        df_all = fish_stim_trains{1,1}(:, stages_ptz(kwid,:));
        corrMat = corrcoef(df_all');
        corrMat(isnan(corrMat)) = 0;
        corrMat = corrMat - diag(diag(corrMat));

        % Proportional Thresholding
        corrMat_thr = threshold_proportional(corrMat, prop_threshold);
        corrMat_bin = double(corrMat_thr > 0);
        
        % Adjacency Matrix for BCT (ensure orientation)
        A = corrMat_bin'; 

        %% Graph Theory Metrics Calculation
        % Global Metrics
        gt.density(fi, kwid)      = density_und(A);
        gt.assortativity(fi, kwid) = assortativity_bin(A, 0);
        gt.modularity(fi, kwid)    = getmodularity(A, swsize);
        gt.efficiency(fi, kwid)    = getefficiency(A, swsize);

        % Map-based Metrics (Node-level)
        gt.gtommap{fi, kwid}       = [ROI_centroids, sum(gtom(A, 1))'];
        
        deg = degrees_und(A);
        gt.degreemap{fi, kwid}     = [ROI_centroids, deg'];
        gt.degree(fi, kwid)        = nanmean(deg);

        bc = betweenness_bin(A);
        gt.betweennessmap{fi, kwid} = [ROI_centroids, bc'];
        gt.betweenness(fi, kwid)    = nanmean(bc);

        % Edge Length Calculation
        edge_len = corrMat_bin .* corrMat; % Binary mask * weights
        gt.edgelengthmap{fi, kwid} = [ROI_centroids, nanmean(edge_len, 2)];
        gt.edgelength(fi, kwid)    = nanmean(gt.edgelengthmap{fi, kwid}(:, 4));

        fprintf('Done.\n');
    end
    
    % Intermittent save for safety
    save('whole_GT24.mat', 'gt', 'N', 'allcond', 'stages_ptz', '-v7.3');
end

%% Helper Functions (Consolidated)
function val = ifthen(cond, a, b)
    if cond, val = a; else, val = b; end
end

function modularity = getmodularity(AdjMat, swsize)
    Ci = zeros(10, swsize);
    for i = 1:10
        rand_idx = randperm(size(AdjMat, 1), swsize);
        Ci(i,:) = modularity_und(AdjMat(rand_idx, rand_idx));
    end
    modularity = mean(nanmedian(Ci, 2));
end

function efficiency = getefficiency(AdjMat, swsize)
    E = zeros(1, 10);
    for i = 1:10
        rand_idx = randperm(size(AdjMat, 1), swsize);
        E(i) = efficiency_bin(AdjMat(rand_idx, rand_idx), 0);
    end
    efficiency = nanmean(E);
end