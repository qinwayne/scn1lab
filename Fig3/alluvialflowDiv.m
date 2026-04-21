function h = alluvialflowDiv(data, cdata,left_labels, right_labels, chart_title,plotsig)
% Copyright 2018 The MathWorks, Inc. 
%
% Plot an alluvial flow diagram.
% left_labels:  Names of categories to flow from.
% right_labels: Names of categories to flow to.
% data:         Matrix with size numel(left_labels) rows by
%               numel(right_labels) columns.
%
% Ideas for future work:
% 1. Get data from a MATLAB table, use table variable names and named rows
%    as labels.
% 2. Interface similar to the plot function, with name-value pairs, optional
%    parameters etc.
    % h = gcf;
    % clf
    % set(h, 'WindowStyle', 'Docked'); % DBG this helps reuse desktop space
    
    % Find axis dimensions and set them
    data_sum = sum(data(:));
    total_gap = 200 - data_sum;
    left_gap_size = total_gap / (size(data, 1)-1);
    right_gap_size = total_gap / (size(data, 2)-1);
    y_height = data_sum + total_gap;
    x_left = 0;
    x_right = 1;
    axis([x_left, x_right, 0, y_height]) % Set limits
    axis ij % origin is top left
    axis off

    patch_colors = [[0.0000, 0.5451, 0.5451];	[0.8549, 0.6471, 0.1255]; 	[0.3333, 0.4196, 0.1843]; [0.6000, 0.1961, 0.8000]];
    
    % grid minor % DBG
    
    hold on
    patch([0 0 1 1], [0 y_height y_height 0], 'w');
    

    % significant_bar
    if plotsig
        plot([0.5 0.5],[0 20],'k','LineWidth',2)
    end
    
    % Plot left categories - one per row
    left_category_sizes = sum(data, 2)';
    
    % These are the top points for each left category, 
    % with gaps added.
    left_category_points = [0 cumsum(left_category_sizes)] + ...
        (0:numel(left_category_sizes)) .* left_gap_size;
    left_category_points(end) = [];
    
    % plot left category bars
    for i = 1:size(patch_colors,1)
        plot([0 0], [left_category_points(i); (left_category_points(i) + left_category_sizes(i))], 'color',patch_colors(i,:), 'LineWidth',4);
    end

    % DBG plot left ticks
    %left_category_tick_starts = zeros(size(left_category_points)) - 0.01;
    %left_category_tick_ends = left_category_tick_starts + 0.02;
    %plot([left_category_tick_starts; left_category_tick_ends], ...
    %     [left_category_points; left_category_points], 'b-');
    
    % Plot right categories - one per column
    right_category_sizes = sum(data, 1);

    % These are the top points for each right category, 
    % with gaps added.
    right_category_points = [0 cumsum(right_category_sizes)] + ...
        (0:numel(right_category_sizes)) .* right_gap_size;
    right_category_points(end) = [];
    
    % plot right category bars
    for i = 1:size(patch_colors,1)
        plot([1 1], [right_category_points(i); (right_category_points(i) + right_category_sizes(i))], 'color',patch_colors(i,:), 'LineWidth',4);
    end

    % DBG plot right ticks
    %right_category_tick_ends = ones(size(right_category_points)) + 0.01;
    %right_category_tick_starts = right_category_tick_ends - 0.02;
    %plot([right_category_tick_starts; right_category_tick_ends], ...
    %     [right_category_points; right_category_points], 'b-');
     

    %
    % Draw the patches, an entire left category at a time
    %
    right_columns_so_far = right_category_points(1:end); % Start at the beginning of each right category and stack as we go.
    patches_per_left_category = size(data, 2);
    normA = abs(cdata) - min(abs(cdata));
    cdatanorm = normA ./ max(normA(:)); % *
    for k_left = 1:size(data, 1) % for each row   
        for k_right = 1:size(data, 2) % for each row
            if cdata(k_left,k_right)>0
                color = 'b';
            else
                color = 'r';
            end
         
        %
        % Calculate the coordinates for all the patches split by the
        % current left category
        %
         
        % Split the left category
        left_patch_points = [0 cumsum(data(k_left, k_right))] + left_category_points(k_left);
        patch_top_lefts = left_patch_points(1:end-1);
        patch_bottom_lefts = left_patch_points(2:end);
        left_category_points(k_left) = patch_bottom_lefts;
         
        % Compute and stack up slice of each right category
        patch_top_rights = right_columns_so_far(1, k_right);
        patch_bottom_rights = patch_top_rights + data(k_left, k_right);
        right_columns_so_far(1,k_right) = patch_bottom_rights;
         
        %
        % Plot the patches
        %
        
        % X coordinates of patch corners
        [bottom_curves_x, bottom_curves_y] = get_curves(0.1, patch_bottom_lefts, 0.9, patch_bottom_rights);
        [top_curves_x,    top_curves_y]    = get_curves(0.9, patch_top_rights,   0.1, patch_top_lefts);        
        X = [ ...
            0;
            % repmat([0; 0], 1, 1); % Top left, bottom left
            bottom_curves_x;
            1;
            % repmat([1; 1], 1, 1); % Bottom right, top right
            top_curves_x
            ];
        
         
        % Y coordinates of patch corners
        Y = [ ...

            (patch_top_lefts+patch_bottom_lefts)/2;
            % patch_top_lefts; 
            % patch_bottom_lefts; 
            bottom_curves_y;
            (patch_bottom_rights+patch_top_rights)/2;
            top_curves_y
            ];

        % if abs(cdata(k_left,k_right))>20
        %     alphaVal = 0.8;
        % else
        %     alphaVal = 0.1;
        % end
        alphaVal = abs(cdatanorm(k_left,k_right));

        patch('XData', X, 'YData', Y, 'FaceColor', color, 'FaceAlpha', alphaVal, 'EdgeColor', 'none');
        end

        % Place left labels
        text(0 - 0.01, ...
            left_category_points(k_left) - left_category_sizes(k_left)./2, ...
            left_labels{k_left},'Color',patch_colors(k_left,:), 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90);

        % Place right labels
        text(1 + 0.01, ...
            right_category_points(k_left) + right_category_sizes(k_left)./2, ...
            right_labels{k_left},'Color',patch_colors(k_left,:), 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Rotation', 90);

    end % for each row

    title(chart_title);
end % alluvialflow

function [x, y] = get_curves(x1, y1, x2, y2)
% x1, x2: scalar x coordinates of line start, end
% y1, y2: vectors of y coordinates of line start/ends
    Npoints = 15;
    t = linspace(0, pi, Npoints);
    c = (1-cos(t))./2; % Normalized curve
    
    Ncurves = numel(y1);
	% Starting R2016b, the following line could be written simply as:
    %   y = y1 + (y2 - y1) .* c';
    y = repmat(y1, Npoints, 1) + repmat(y2 - y1, Npoints,1) .* repmat(c', 1, Ncurves);
    x = repmat(linspace(x1, x2, Npoints)', 1, Ncurves);
end  % get_curve