function [curr] = scatterfit(var1, var2, legendAdd, varargin)
% var1 should be ndataPoints x nSubjects or conditions,etc.
% var2 should have same structure

% uses a number of default colours, if nConditions exceeds that then it
% uses the last one 

cols = {[230/255, 159/255, 1/255],   [86/255, 180/255, 233/255], ...
        [1/255,   158/255, 115/255], [240/255, 228/255, 66/255], ...
        [1/255,   114/255, 178/255], [213/255, 94/255, 1/255], ...
        [204/255, 121/255, 167/255]};

nCols = size(var1, 2);

grayVec  = [0.5804, 0.5843, 0.600];
greenVec = [0.3294    0.6314    0.5490];
%greenVec = grayVec;

cols{1} = greenVec;

corrs = [];
for i = 1:nCols
    
    if i > size(cols, 2)
        cols = cols{end};
    end
    
    % first remove nans
    X = var1(:, i); X = var1(~isnan(var1(:, i)), i);
    Y = var2(:, i); Y = var2(~isnan(var2(:, i)), i);
    
    if ~isempty(varargin)
        if length(varargin) == 2
            [a,b] = corr(X, Y, 'type', 'spearman');
        end
    else
        [a,b] = corr(X, Y);
    end
    
    if ~isempty(varargin)       
        fo(i)=scatter(X, Y, 70, 'filled', 'MarkerEdgeColor', [1/255 1/255 1/255],...
        'MarkerFaceColor', varargin{1},...
        'LineWidth', 1.5); hold on
    else
        fo(i)=scatter(X, Y, 70, 'filled', 'MarkerEdgeColor', [1/255 1/255 1/255],...
        'MarkerFaceColor', cols{i},...
        'LineWidth', 1.5); hold on
    end

    h = lsline;

    if ~isempty(varargin)       
        set(h(1),'color',varargin{1}, 'linewidth', 1.5)
    else
        set(h(1),'color',cols{1}, 'linewidth', 1.5)
    end

    corrs(i, 1:2) = [a, b];
end

% keyboard

leg = cell(nCols, 1);
for i = 1:nCols
    leg{i} = strcat(sprintf('r = %.2f (%.2f)', corrs(i, :)));
end

% between each add an empty space 
if legendAdd == 1
    legend(fo, leg, 'location', 'best'); 
end

%if ~isempty(varargin)    
%    figElements(curr, varargin{:});
%end

box off

end
