%% SCN1A Density Comparison Analysis
% This script compares WT and MUT fish density data across 3D pixel boxes,
% performs statistical testing (Mann-Whitney U), computes pixelwise differences,
% and generates a video visualization of the results.

%% Setup
clear all;
cd("D:\OneDrive - The University of Melbourne\Research\SCN1\nodedensity")
load("density_comparison2.mat");  % Load preprocessed density data

finalMat = zeros(targetSize);     % Initialize output matrix
pixel_add = 1;                    % Padding around each voxel
all_fish = [WT_fish SCN1_fish];   % Combine WT and MUT fish indices

%% Statistical Comparison: Mann-Whitney U Test
% For each voxel, test whether WT and MUT densities differ significantly

for xi = 1+pixel_add : targetSize(1)-pixel_add
    for yi = 1+pixel_add : targetSize(2)-pixel_add
        for zi = 1+pixel_add : targetSize(3)-pixel_add

            wt_mov_den = [];  % WT densities
            mut_mov_den = []; % MUT densities

            % Collect WT voxel data
            for fi = 1:length(WT_fish)
                data = reszden{fi}(xi-pixel_add:xi+pixel_add, ...
                                  yi-pixel_add:yi+pixel_add, ...
                                  zi-pixel_add:zi+pixel_add);
                wt_mov_den = [wt_mov_den; data(:)];
            end

            % Collect MUT voxel data
            for fi = length(WT_fish)+1 : length(all_fish)
                data = reszden{fi}(xi-pixel_add:xi+pixel_add, ...
                                  yi-pixel_add:yi+pixel_add, ...
                                  zi-pixel_add:zi+pixel_add);
                mut_mov_den = [mut_mov_den; data(:)];
            end

            % Perform Mann-Whitney U test
            try
                [p, h] = ranksum(wt_mov_den, mut_mov_den);
                if p < 0.05
                    finalMat(xi, yi, zi) = nanmean(wt_mov_den) - nanmean(mut_mov_den);
                end
            end
        end
    end
end

%% Pixelwise Mean Subtraction Between KO and WT
load('../Warp analysis/warp_data.mat');  % Load warped density data

zrange = [1, size(KO_data,3)];
yrange = [1, size(KO_data,2)];
xrange = [1, size(KO_data,1)];

for plane = zrange(1):zrange(2)
    for x = xrange(1):xrange(2)
        for y = yrange(1):yrange(2)
            KO_temp = nanmean(double(KO_data(x,y,plane,:)));
            WT_temp = nanmean(double(WT_data(x,y,plane,:)));
            pixelwise_diff(x,y,plane) = KO_temp - WT_temp;
        end
    end
end

%% Video Visualization of Slice-by-Slice Comparison
close all;
v = VideoWriter('Video S1.avi');
v.FrameRate = 10;
open(v);

figure('Position', [100, 100, 1600, 1000]);
uplim = max(pixelwise_diff,[],'all');

scalex = (850-0)/targetSize(1);
scaley = (430-0)/targetSize(2);
oldcanvas = [];

for i = 1:size(pixelwise_diff,3)
    if mod(i, 14) == 10 && i < 248
        % Left subplot: statistical comparison heatmap
        subplot(1,2,1);
        if ~isempty(oldcanvas)
            hold off;
            imagesc(oldcanvas, 'AlphaData', canvasAlpha);
            hold on;
        end

        data_temp = squeeze(rot90(finalMat(:,:, (i+4)/14), 2));
        imAlpha = ones(size(data_temp));
        imAlpha(data_temp == 0) = 0;  % Transparent where no difference

        % Create white canvas and center the image
        canvas = ones([150, 120]);
        startY = floor((150 - size(data_temp, 1)) / 2);
        startX = floor((120 - size(data_temp, 2)) / 2) + 1;
        canvas(startY:startY+size(data_temp,1)-1, ...
               startX:startX+size(data_temp,2)-1) = data_temp;

        canvasAlpha = zeros(size(canvas));
        canvasAlpha(startY:startY+size(imAlpha,1)-1, ...
                    startX:startX+size(imAlpha,2)-1) = imAlpha;

        imagesc(canvas, 'AlphaData', canvasAlpha);
        oldcanvas = canvas;

        clim([-0.3 0.3]);
        axis off;
        colormap(flip(redblueRealTecplot));
        set(gca, 'box','off', 'TickDir','out', ...
            'XTickLabel',[], 'XTick',[], ...
            'YTickLabel',[], 'YTick',[], ...
            'xcolor',[1 1 1], 'ycolor',[1 1 1]);
    end

    % Right subplot: KO vs WT pixelwise difference
    subplot(1,2,2);
    imagesc(rot90(pixelwise_diff(:,:,i),1), [-uplim, uplim]);
    colormap(redblueRealTecplot);
    axis off;
    xlim([0, size(pixelwise_diff,1)+10]);
    ylim([0, size(pixelwise_diff,2)]);
    set(gca, 'box','off', 'TickDir','out', ...
        'XTickLabel',[], 'XTick',[], ...
        'YTickLabel',[], 'YTick',[], ...
        'xcolor',[1 1 1], 'ycolor',[1 1 1]);

    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);

%% Helper Function: Compute ROI Density in 3D Cubes
function [N, avgMat] = getROIsDensity(data, cubesize)
% Computes the number of points in each 3D cube of size `cubesize`
% Returns:
%   N       - number of data points
%   avgMat  - 3D matrix of counts per cube

N = length(data);
xVal = min(data(:,1)); xi = 1;

while xVal < max(data(:,1))
    xidx = data(:,1) >= xVal & data(:,1) < xVal + cubesize;
    yVal = min(data(:,2)); yi = 1;

    while yVal < max(data(:,2))
        yidx = data(:,2) >= yVal & data(:,2) < yVal + cubesize;
        zVal = min(data(:,3)); zi = 1;

        while zVal < max(data(:,3))
            zidx = data(:,3) >= zVal & data(:,3) < zVal + cubesize;
            roiidx = find(xidx & yidx & zidx);

            if ~isempty(roiidx)
                avgMat(xi, yi, zi) = length(roiidx);
            end

            zVal = zVal + cubesize;
            zi = zi + 1;
        end

        yVal = yVal + cubesize;
        yi = yi + 1;
    end

    xVal = xVal + cubesize;
    xi = xi + 1;
end
end
