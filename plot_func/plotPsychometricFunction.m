function curr = plotPsychometricFunction(inp, curr, figInfo, varargin)
% inp takes a struct from the plotPsyche package 

greenVec = [0.3294, 0.6314, 0.5490];
blueVec   = [1/255,   114/255, 178/255];

    
    

plot([0, 0], [0, 1], ':', 'linewidth', 0.5, 'color', [0.05, 0.05, 0.05]); hold on
plot([min(inp.x), max(inp.x)], [0.5 0.5], ':', 'linewidth', 0.5, 'color', [0.05, 0.05, 0.05])


if ~isempty(varargin)
    colVec = varargin{1};
    plot(inp.curve(:, 1), inp.curve(:, 2), 'color', colVec, 'linewidth', 3); hold on
else
    plot(inp.curve(:, 1), inp.curve(:, 2), 'color', blueVec, 'linewidth', 3); hold on
end

box off 
scatter(inp.x, inp.y, 'filled', 'MarkerEdgeColor', 'k', ...
                                'MarkerFaceColor', greenVec, ...
                                'LineWidth', 1.5)
                               
figElements(curr, figInfo{:})


end