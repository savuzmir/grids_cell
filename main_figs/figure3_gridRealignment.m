function figure3_gridRealignment()
clear;clf;

%% load data structs
% cd(go into data folder) 
rel_files = {'auxStruct', 'out_grid_realignment', 'out_grid_realignment_distance', 'out_grid_realignment_distance_control'};
for fl = 1:length(rel_files)
    load(rel_files{fl});
end

%% prepare and plot data
% deflate containers
withinStimBetas  = removenans(out.withinStimBetas, 1);
betweenStimBetas = removenans(out.betweenStimBetas, 1);
withinStimBetasOnoff  = removenans(out.withinStimBetasOnoff, 1);
betweenStimBetasOnoff = removenans(out.betweenStimBetasOnoff, 1);

curr = figure(11);

subplot(3, 2, 1)
figInfo = {sprintf('VMPFC_{StimSet A -> StimSet A}'), 'Symmetries', 'Betas [a.u.]', [], [-0.03, 0.15], [1 2 3 4 5], ...
    {'4', '5', '6', '7', '8'}, [], {}, 24, [], []};
[curr] = plotBar(withinStimBetas', figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, [0, 1], 0, 0);

[a,b,c,d]=ttest(withinStimBetas')
% [0.1292 -1.5157 4.7673 1.1682 -1.1998]
% 0.9004    0.1681    0.0014    0.2764    0.2645
out.statsWithin = [a; b];
out.withinBetas = withinStimBetas;
subplot(3, 2, 2)
figInfo = {sprintf('VMPFC_{StimSet A -> StimSet B}'), 'Symmetries', '', [], [-0.03, 0.15], [1 2 3 4 5], ...
    {'4', '5', '6', '7', '8'}, [], {}, 24, [], []};
[curr] = plotBar(betweenStimBetas', figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, [0, 1], 0, 0);

[a, b, c, d]=ttest(betweenStimBetas')
% [-2.7434 -0.4079 -0.0576 -1.6424 1.3475]
% 0.0151    0.6891    0.9548    0.1213    0.1978
out.statsBetween = [a; b];
out.betweenBetas = betweenStimBetas;

subplot(3, 2, 3)
figInfo = {sprintf(''), 'Angles', 'z-score [a.u.]', [], [-0.2, 0.2], [1 3 5 7 9 11], ...
    {'\omega', '\omega + 60', '\omega + 120', '\omega + 180', '\omega + 240', '\omega + 300'}, [], {}, 24, [], []};

[curr] = plotBar(withinStimBetasOnoff', figInfo, auxStruct.pinkVec, auxStruct.grayVec, curr, [0, 1], 0, 0);

figInfo = {sprintf(''), 'Angles', '', [], [-0.2, 0.2], [1 3 5 7 9 11], ...
    {'\omega', '\omega + 60', '\omega + 120', '\omega + 180', '\omega + 240', '\omega + 300'}, [], {}, 24, [], []};
subplot(3, 2, 4)
[curr] = plotBar(betweenStimBetasOnoff', figInfo, auxStruct.pinkVec, auxStruct.grayVec, curr, [0, 1], 0, 0);

trueDiff    = statsOrientation.permDiff(1, 2, 1) - statsOrientation.permDiff(1, 1, 1);
permAvgDiff = statsOrientation.permDiff(1, 2, 2:end) - statsOrientation.permDiff(1, 1, 2:end);

subplot(3, 2, 6)
perm_histogram(squeeze(permAvgDiff), {repmat(trueDiff, [round(length(permAvgDiff)/7),  1]), 1:round((length(permAvgDiff)/7))}, {'edgecolor', 'none', 'facecolor', auxStruct.grayVec}, {'LineWidth', 3', 'color', 'k'}); hold on
figElements(curr, 'More similar within stimuli sets', '(Between - Within) grid orientation difference', 'Frequency', [], [], [], ...
    {}, [], {}, 24, [], []);
box off
legend('Null distribution', '(Between - Within) pairwise distance', 'Location', 'northoutside'); legend box off

%% supplements; control analysis

curr = figure(7)

trueDiff    = statsOrientation_control.permDiff(1, 2, 1) - statsOrientation_control.permDiff(1, 1, 1);
permAvgDiff = statsOrientation_control.permDiff(1, 2, 2:end) - statsOrientation_control.permDiff(1, 1, 2:end);

perm_histogram(squeeze(permAvgDiff), {repmat(trueDiff, [round(length(permAvgDiff)/7),  1]), 1:round((length(permAvgDiff)/7))}, {'edgecolor', 'none', 'facecolor', auxStruct.grayVec}, {'LineWidth', 3', 'color', 'k'}); hold on
figElements(curr, 'More similar within stimuli sets', '(Between - Within) grid orientation difference', 'Frequency', [], [], [], ...
    {}, [], {}, 24, [], []);
box off
legend('Null distribution', '(Between - Within) pairwise distance', 'Location', 'northoutside'); legend box off


end