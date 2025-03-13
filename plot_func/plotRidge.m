function [curr, minVals] = plotRidge(input, overlap, plotType, cols, curr)
% INPUTS:
% INPUT (MATRIX): state by time by features (features are either channels or units)
% OVERLAP(FLOAT): determines degree of overlap between individual plots
% PLOTTYPE(INT): yields a prob density-like (0) or plotmse like ridge (1) plot
% COLS (MATRIX): num states by col matrix

% OUTPUT:
% curr - figure handle
% minVals - values used to adjust y axis ticks

% Sample points
N = size(input, 2);
nStates = size(input, 1);
% Plot options
mini = 1;
maxi = N;
overl = overlap;

% data parse
avgData = nanmean(input, 3)';

semData = sem(input, 3)';

lb_se = avgData - semData;
ub_se = avgData + semData;

% Get the position of each dataset
y = cumsum(max(avgData,[],1))*(1-overlap);


xLen = linspace(mini, maxi, N);
minVals = [];

for ii = nStates:-1:1
    minVal = min(avgData(:, ii)) + y(ii);
    
    minVals(ii) = minVal;
    
    if plotType == 1
    
        filling = fill([xLen fliplr(xLen)], [ub_se(:,ii)' + y(ii) fliplr(lb_se(:, ii)' + y(ii))], cols(ii,:), 'linestyle', 'none'); hold on
    else    
        filling = fill([xLen fliplr(xLen)], [avgData(:, ii)' + y(ii), fliplr(repmat(minVal, [1, length(xLen)]))], cols(ii,:), 'linestyle', 'none'); hold on
    end
    
    plot(xLen, repmat(minVal, [length(xLen), 1]), 'k');
    plot(xLen, avgData(:,ii)+ y(ii), 'k','LineWidth',1); hold on
end
