function newStats = figure_GridOrientation_v3_supplements(chanInfo, all_units, dayInformationStimset, nPerm, auxStruct, reg)

reg = 4;
grayVec = auxStruct.grayVec;
symmetries = 6; % this will select only the sixfold symmetry for testing
control    = 0;
[statsOrientation] = computingOrientationDistance_pap(chanInfo, all_units, dayInformationStimset, nPerm, symmetries, auxStruct, reg, control);

symmetries = 4:8;
[curr, out] = predictingAcrossRegionsBetas_pap(chanInfo, all_units, auxStruct, symmetries, dayInformationStimset, reg);

trueDiff    = statsOrientation.permDiff(1, 2, 1) - statsOrientation.permDiff(1, 1, 1);
permAvgDiff = statsOrientation.permDiff(1, 2, 2:end) - statsOrientation.permDiff(1, 1, 2:end);

subplot(3, 2, 6)
perm_histogram(squeeze(permAvgDiff), {repmat(trueDiff, [round(length(permAvgDiff)/7),  1]), 1:round((length(permAvgDiff)/7))}, {'edgecolor', 'none', 'facecolor', grayVec}, {'LineWidth', 3', 'color', 'k'}); hold on
figElements(curr, 'More similar within stimuli sets', '(Between - Within) grid orientation difference', 'Frequency', [], [], [], ...
    {}, [], {}, 24, [], []);
box off
legend('Null distribution', '(Between - Within) pairwise distance', 'Location', 'northoutside'); legend box off

confInt = [0, 0.975];
distances = nan(120, 2);

withinVec   = statsOrientation.withinVecFull;
betweenVec  = statsOrientation.betweenVecFull;

distances(1:length(withinVec), 1) = withinVec;
distances(1:length(betweenVec), 2) = betweenVec;

% [a,b,~, stats] = ttest2(withinVec, betweenVec);

% subplot(3, 2, 1)
% withinCol  = [0.3294, 0.5314, 0.5490];
% betweenCol = [213/255, 94/255, 1/255];
% 
% figInfo = {'Average pairwise grid orientation distance', 'Distance', 'Distance [deg]', [], [0 30], [1 2], ...
%            {'Within', 'Between'}, [], {}, 28, [], []};
% [curr] = plotBar(distances, figInfo, withinCol, betweenCol, curr, confInt, 0);

% [a,b,c,d]=ttest2(withinVec, betweenVec)
% t = 2.4, p = 0.02

currFields = fields(statsOrientation);

% rewrite stats into final stat container
for fl = 1:length(currFields)
    newStats(reg).(currFields{fl}) = statsOrientation.(currFields{fl});
end

newStats(reg).statsWithin  = out.statsWithin; 
newStats(reg).statsBetween = out.statsBetween;

set(curr, 'Units', 'normalized', 'Position', [0, 0, 0.2, .8])
set(curr, 'visible', 'off')
curr.Units = 'pixels';
curr.OuterPosition = [0 0 14500, 12000];
res = 600;
set(curr, 'PaperPositionMode', 'manual');
curr.PaperUnits = 'inches';
curr.PaperPosition = [0, 0, 14500, 12000]/res;

saveImage(curr, ['Fig3_new - VMPFC pairwise difference sess', num2str(reg)], 'C:\Users\sebas\OneDrive - University College London\Sebastijan\KennerleyLab\Figures\grids_working\paper_figures\revision\', [14500, 12000])

clear curr

% run control analysis
grayVec = auxStruct.grayVec;
symmetries = 6; % this will select only the sixfold symmetry for testing
[statsOrientation] = computingOrientationDistance_pap(chanInfo, all_units, dayInformationStimset, nPerm, symmetries, auxStruct, reg);

end