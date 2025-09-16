%% Repeated Measures ANOVA with Bonferroni Correction
% This script loads multiple scenario datasets, normalizes correlation distances,
% performs repeated measures ANOVA across conditions, and visualizes results
% with significance markers and fading effects in a video.

close all; clear; clc;
cd("D:\OneDrive - The University of Melbourne\Research\SCN1\paper\power_law_2018")

%% Load Scenario Files
scenarioFiles = {
    'main20_distance_comp_stat_wtn_reg2reg.mat', ...
    'main20_distance_comp_stat_btn_reg2reg.mat', ...
    'main20_distance_comp_stat_wtn_reg2hem.mat', ...
    'main20_distance_comp_stat_btn_reg2hem.mat'
};

distint = 30;
colors = getColours();

%% Setup Video Writer
v = VideoWriter('Video S2.mp4','MPEG-4');
v.FrameRate = 1;
v.Quality = 100;
open(v);

%% Define Parameters
num_lengths = size(pMat, 1);
num_scenarios = size(pMat, 2);
num_regions = size(pMat{1,1}, 1);
num_timepoints = size(pMat{1,1}, 2);

regionLabels = {'Telencephalon','Diencephalon','Mesencephalon','Rhombencephalon'};
scenarioLabels = {'Ipsi: Reg→Reg','Cont: Reg→Reg','Ipsi: Hem→Reg','Cont: Hem→Reg'};
regioncols = getColours();

% Generate distance bin labels
for di = 2:size(corrdisMat,1)
    dislabels{di-1} = mat2str(distint*(di-1));
end

% Generate time labels
kwid = 1;
winsize = 300;
slidewin = winsize / 2;
while 700 + slidewin*(kwid-2) + winsize < 4350
    starttime = (kwid-1)*slidewin + 100;
    endtime = kwid*slidewin + winsize + 100;

    if starttime < 710
        timelabels{kwid} = 'Pre PTZ';
    else
        timelabels{kwid} = [mat2str(starttime/2 - 300) 's'];
    end
    kwid = kwid + 1;
end

%% Generate Visualization Frames
fadeDuration = 4;
fadeFrames = fadeDuration * v.FrameRate;

for kwid = 1:num_timepoints
    for f = 1:fadeFrames
        fadeAlpha = 1;  % You can change to f / fadeFrames for true fade-in

        figure(1);
        set(gcf, 'Position', [100, 100, 2400, 1100]);

        sgtitle(['Time Step: ' timelabels{4+kwid}], 'FontSize', 16, 'FontWeight', 'bold');

        % Progress bar
        progress = kwid / num_timepoints;
        barHeight = 0.01;
        barY = 0.01;
        barColor = [0.2+kwid*0.6/num_timepoints, 0.6, 0.8-kwid*0.6/num_timepoints];

        annotation('textbox', [0.03, barY, 0.15, 0.1], ...
            'String', 'PTZ added', 'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'FontWeight', 'bold', 'Color', [0, 0, 0]);

        annotation('arrow', [0.06 0.06], [0.1 0.03], 'Color', [0, 0, 0], 'LineWidth', 1.5);
        annotation('rectangle', [0, barY, 1, barHeight], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
        annotation('rectangle', [0, barY, progress, barHeight], 'FaceColor', barColor, 'EdgeColor', 'none');

        % Plot each region
        for regi = 1:num_regions
            subplot(2,2,regi);
            if f > 1, cla; end
            hold on;

            for s = 1:num_scenarios
                yOffset = 5 * (num_scenarios - s + 1);

                for i = 1:num_lengths
                    val = meanMat{i,s}(regi,kwid);
                    pval = pMat{i,s}(regi,kwid);

                    xPos = i * 2 - 1;
                    yPos = yOffset;

                    colors = flipud(redblueRealTecplot(100));
                    val_clipped = max(min(val, 1), -1);
                    idx = round((val_clipped + 1)/2 * (size(colors,1)-1)) + 1;
                    rgb = colors(idx,:);

                    scatter(xPos, yPos, -log(pval)*60, rgb, 'filled', 'MarkerFaceAlpha', fadeAlpha);

                    if ~isempty(txtMat{i,s}{regi,kwid}) && f < fadeFrames
                        scatter(xPos, yPos, -log(pval)*60, [0 0 0], ...
                            'MarkerEdgeAlpha', fadeAlpha, 'LineWidth', 2);
                    end

                    text(-1, yOffset, scenarioLabels{s}, 'FontSize', 10, ...
                        'HorizontalAlignment', 'right');
                end

                title(regionLabels{regi}, 'FontSize', 14, 'Color', regioncols{regi});
                xlabel('Distance bin (\mum)');
                xlim([0 35]);
                xticks(1:2:35);
                xticklabels(dislabels);
                ylim([0 25]);
                yticks([]);
                grid on;
            end
        end

        drawnow;
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
end

close(v);
