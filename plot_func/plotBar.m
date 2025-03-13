function [curr] = plotBar(input, figInfo, colorVec1, colorVec2, curr, confInt, indivPoints, paired_samples, varargin)
% [curr] = plotBar(input, figInfo, colorVec1, colorVec2, curr, confInt, indivPoints, varargin)

full_cols = {[230, 159, 1],  [86, 180, 233], [1, 158, 115], ...
             [240, 225, 66], [1, 114, 178],  [213, 94, 1], ...
             [141, 211,199], [204, 121, 167],[148, 149, 153], ...
             [84, 161, 140], [255, 255, 179],[190, 186, 218], ...
             [251, 128, 114],[128, 177, 211],[253, 180, 98], ...
             [179, 222, 105],[252, 205, 229],[217, 217, 217], ...
             [188, 128, 189],[204, 235, 197]};

full_cols = cellfun(@color_norm, full_cols, 'UniformOutput', false);

palette = [0.9, 0.3, 0.3; 
           0.9, 0.6, 0.1; 
           0.9, 0.9, 0.1; 
           0.1, 0.9, 0.1; 
           0.4, 0.4, 0.9; 
           0.4, 0.1, 0.5;
           0.9, 0.2, 0.7;
           0.7, 0.3, 0.3;
           0.1, 0.1, 0.1;
           0.5, 0.5, 0.5;
           0.7, 0.7, 0.9;
           0.7, 0.9, 0.7];

palette = viridis;
palette = palette(1:6:end, :);

full_cols = palette;

rel_col = [colorVec1; colorVec2];

if indivPoints

    % get non-nan number of inputs
    if paired_samples
        tmp = input(:, 1);
        tmp = tmp(~isnan(tmp));
        nStates = size(tmp);
    else
        tmp = input(1, :);
        tmp = tmp(~isnan(tmp));
        nStates = size(tmp, 2);
    end

    for i = 1:nStates

        if paired_samples

            vals = [];
            X = find(~isnan(input(i, :)));
            
            for e = 1:length(X)
                vals(e) = rand * 0.15 + 0.025;
            end
            X = X + vals;
            Y = input(i, :);
            Y = Y(~isnan(Y));
        else
            X = zeros(sum(~isnan(input(:, i))), 1) + i;
            Y = input(:, i);
            Y = Y(~isnan(Y));
            X = X';
            Y = Y';
            
            vals = [];
            % add noise
            for e = 1:length(X)
                vals(e) = rand * 0.15 + 0.025;
            end            
            X = X + vals;

        end

        alph_lev = 0.7;

%         if colored
%             colo = [full_cols(i, :)];
%         else
%             colo = 'k';
%         end
        
        if isempty(Y)
            continue
        end

        if paired_samples
            line(X, Y, 'Color', [0, 0, 0, 0.3], 'LineWidth', 1); hold on
        end

        for e = 1:size(X, 2)
            if size(input, 2) > 2
                rel_curr_col = rel_col(1, :);
            else
                rel_curr_col = rel_col(e, :);
            end

            scatter(X(e), Y(e), 'filled', ... 
                      'MarkerEdgeColor', [0.2 0.2 0.2],     ...
                      'MarkerFaceColor', rel_curr_col, ...
                      'MarkerFaceAlpha',  alph_lev, ...
                      'MarkerEdgeAlpha', alph_lev, ... 
                      'LineWidth', 1); hold on
        end

    end
end

nStates = size(input, 2);
b = bar(nanmean(input, 1), 'edgeColor', 'k', 'LineWidth', 3); hold on
b.FaceColor = 'flat';

if size(colorVec1, 2) == 4
    colAlpha = colorVec1(1, 4);
    colorVec1 = colorVec1(:, 1:3);
    colorVec2 = colorVec2(:, 1:3);
    
else
    colAlpha = 0.7;
end

b.FaceAlpha = colAlpha;

if all(colorVec1==colorVec2) & size(input, 2) > 1 % means you don't want to color them separately so we can use nStates
    
    b.CData(:, :) = repmat(colorVec1, [nStates, 1]);
    
elseif any(colorVec1~=colorVec2) && isempty(varargin) & size(input, 2) > 1 % means you want to color st 1 as col 1, st 2 as col 2 st 3 as col 1 etc.
%     keyboard

try
    b.CData(1:2:end, :) = repmat(colorVec1, [nStates/2, 1]);
    b.CData(2:2:end, :) = repmat(colorVec2, [nStates/2, 1]);
catch
    b.CData(1, :) = repmat(colorVec1, [1, 1]);
    b.CData(2, :) = repmat(colorVec2, [1, 1]);
end

elseif size(input, 2) == 1
    b.CData(1, :) = repmat(colorVec1, [1, 1]);
else
    
    % we select which one we want to color
    b.CData(1:nStates, :) = repmat(colorVec2, [nStates, 1]);
    b.CData(varargin{1}, :) = repmat(colorVec1, [1, 1]);
    
end

if confInt(1)
    w = norminv(confInt(2));
    err = w*nanstd(input, [], 1)./sqrt(sum(~isnan(input), 1));
else
    err = nanstd(input, [], 1)./sqrt(sum(~isnan(input), 1));
end

errorbar(1:nStates, nanmean(input, 1), err,  '.', 'CapSize', 0, 'linewidth', 3, 'color', [1/255, 1/255, 1/255])

figElements(curr, figInfo{:})
box off
end
