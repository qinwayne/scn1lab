%% Initial Setup and Parameters
clc; clear; close all;

% Group Definitions
WT_fish = {'07','14','17','19','20','24','25','27','32','35','42','58','59','37','52','53','63'};
SCN1_fish = {'01','05','08','16','21','22','38','46','47','45','60','61','62','64','65','66','67','68'};
all_fish = [WT_fish, SCN1_fish];

% Constants
fs = 2;             % Imaging sampling frequency
fsv = 199.94;       % Tail video sampling frequency
min_sz_dur = 1*fs;  % Minimum seizure duration (frames)
winsz = 40;
slidsz = winsz/2;
tailseg = 3:5;      % Tail segments to analyze

% Paths
base_data_path = '..\data\SCN1LabData\';
raw_data_path  = 'D:\Wei\SCN1LabData\';

% Pre-allocate storage
num_fish = length(all_fish);
motion_corr     = cell(num_fish, 1);
motor_regressor = cell(num_fish, 1);
seizureTime_mc  = cell(num_fish, 1);
seizureTime_ci  = cell(num_fish, 1);
seizureTime_tm  = cell(num_fish, 1);
seizureSkewness = cell(num_fish, 1);

%% Generate Sliding Windows
i = 1;
while 1 + slidsz*(i-1) + winsz < 4200
    stages_ptz(i,:) = 1+slidsz*(i-1) : 1+slidsz*(i-1)+winsz;
    allcond{i} = sprintf('stage%d', i-1);
    i = i + 1;
end
allcond{1} = 'pre';

%% Main Processing Loop
for fi = 1:num_fish
    curr_fish = all_fish{fi};
    fprintf('Processing Fish %d/%d: %s\n', fi, num_fish, curr_fish);
    
    % 1. Load Suite2p Motion Data
    fish_dir = dir(fullfile(base_data_path, 'suite2p', ['suite2p_*fish' curr_fish '_*']));
    if isempty(fish_dir), continue; end
    
    suite2p_path = fullfile(fish_dir.folder, fish_dir.name, 'plane15', 'Fall.mat');
    load(suite2p_path, 'ops');
    
    motion_corr{fi} = ops.corrXY;
    motor_regressor{fi} = ((-motion_corr{fi}) * 100) + 0.9; 
    motor_detrend = detrend(motor_regressor{fi});

    % Identify peaks in motion correction
    peak_locs_mc = find(isoutlier(diff(motor_detrend), "percentiles", [1 99]));
    seizureTime_mc{fi} = peak_locs_mc([true, diff(peak_locs_mc) > min_sz_dur]) / fs;

    % 2. Load Raw DFF Data
    geno = ifthen(fi > length(WT_fish), 'Hom', 'WT');
    raw_path = fullfile(raw_data_path, geno, ['raw_fish_std_fmt_' curr_fish '.mat']);
    load(raw_path, 'fish_stim_trains');
    
    meandff = nanmean(fish_stim_trains{1,1});
    peak_locs_ci = find(isoutlier(diff(meandff), "percentiles", [1 99]));
    szlocs = peak_locs_ci([true, diff(peak_locs_ci) > min_sz_dur]) / fs;
    seizureTime_ci{fi} = szlocs;

    % 3. Calculate Network Skewness
    ns = NaN(1, length(szlocs));
    for k = 1:length(szlocs)
        if szlocs(k) > 300
            idx = round(szlocs(k) * fs);
            % Bound checking for window
            if idx + 20 > 4200
                df_window = fish_stim_trains{1,1}(:, idx-20 : end);
            else
                df_window = fish_stim_trains{1,1}(:, idx-20 : idx+20);
            end
            
            cMat = corrcoef(df_window');
            cMat(isnan(cMat)) = 0;
            cMat = cMat - diag(diag(cMat)); % Remove self-correlation
            ns(k) = skewness(cMat, 1, 'all');
        end
    end
    seizureSkewness{fi} = ns;

    % 4. Load and Process Tail Movement
    fn_base = fullfile(base_data_path, 'tail', ['Fish_' curr_fish '_base_beh.h5']);
    fn_ptz  = fullfile(base_data_path, 'tail', ['Fish_' curr_fish '_5mM_beh.h5']);
    
    figure('Name', ['Fish ' curr_fish]); hold on;
    t_img = (0:length(meandff)-1)/fs;
    
    if isfile(fn_base) && isfile(fn_ptz)
        extraframe = 6000;
        
        % Read H5 Data
        theta_base = h5read(fn_base, '/tail_theta');
        theta_ptz  = h5read(fn_ptz, '/tail_theta');
        
        data_base = sum(abs(theta_base(2000:end-4001, tailseg)), 2);
        data_ptz  = sum(abs(theta_ptz(1:end-extraframe, tailseg)), 2);
        
        t_base_vec = (0:size(data_base, 1)-1) / fsv;
        t_ptz_vec  = (1:size(data_ptz, 1)) / fsv + t_base_vec(end);
        
        full_tail_data = [data_base; data_ptz];
        full_t_vec = [t_base_vec, t_ptz_vec];
        mean_tail_data = movmean(full_tail_data, 10, 'omitmissing');
        
        % Tail Seizure Detection
        peak_locs_tm = find(isoutlier(mean_tail_data, "percentiles", [0.1 99.9])) / fsv;
        seizureTime_tm{fi} = peak_locs_tm([true; diff(peak_locs_tm(:)) > min_sz_dur/fs]);
        
        % Plot Tail Data
        plot(full_t_vec, zscore(mean_tail_data), 'DisplayName', 'Avg Tail Movement');
        plot(t_ptz_vec, zscore(data_ptz), 'g', 'DisplayName', 'Tail Movement (PTZ)');
        plot(seizureTime_tm{fi}, -10, 'g.', 'HandleVisibility', 'off');
    end

    % 5. Visualization Cleanup
    plot(t_img, 10 + zscore(motor_detrend), 'b', 'DisplayName', 'Motion Corr');
    plot(seizureTime_mc{fi}, 9, 'b*', 'HandleVisibility', 'off');
    plot(t_img, zscore(meandff) + 20, 'r', 'DisplayName', 'dF/F');
    plot(seizureTime_ci{fi}, 19, 'r*', 'HandleVisibility', 'off');
    
    legend('Location', 'best');
    title(['Activity Summary: Fish ' curr_fish]);
end

save('seizure_time.mat', 'seizureTime_ci', 'seizureTime_mc', 'seizureTime_tm', 'seizureSkewness');

%% Helper Functions
function val = ifthen(cond, a, b)
    if cond, val = a; else, val = b; end
end

function newdta = multidetrend(data)
    win = 30;
    meandata = zeros(size(data));
    for i = 1:length(data)
        curr_win = min(win, length(data)-i);
        windata = data(i:i+curr_win);
        % Logic: mean of non-outliers in the detrended window
        meandata(i) = mean(windata(~isoutlier(detrend(windata))));
    end
    newdta = data - meandata;
end