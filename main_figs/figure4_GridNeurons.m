function out = figure4_GridNeurons(auxStruct)
clear; clf;

%% load data structs
%cd(go into data folder) 
rel_files = {'unitPhaseLocked', 'avgTheta', 'auxStruct', 'AuxiliaryCells', 'stored_betas_roi', 'cpd_inf_perm_vmpfc', 'gridPhaseBetas', 'averagedFRInformation', 'averagedPhaseInformation', 'permutationStats_vmpfc_cells', 'sample_chan_overallPhases'};
for fl = 1:length(rel_files)
    load(rel_files{fl});
end

%% ---------------------------
% start cell figure generation
% top row: phase + oscillation concept, cell average across phases, two sample cells

curr = figure(66);

subplot(6, 6, 1)
[dat, f]     = cwt(avgTheta(141, :), auxStruct.Fs, 'amor');

angle_inp    = angle(dat);

% take one half and mirror it
samp_th          = 1073:1146;
p_a              = avgTheta(141, samp_th); % sample chan w theta 
full_theta_cycle = [avgTheta(141, samp_th)  fliplr(avgTheta(141, samp_th(2:end)))]; % start from 2nd index for right hand part of cycle

samp_theta_full_len = 1073:(1072 + length(samp_th)*2 - 1);

plot(normalize_bound(full_theta_cycle, -2, 3.5), 'color', auxStruct.blueVec, 'linewidth', 3); hold on
plot(angle(dat(f>6.61 & f < 7.08, samp_theta_full_len)), 'color', 'k', 'linewidth', 2.5); box off % we want to capture 1 cycle defined within this freq
figElements(curr, '', '', '', [], [], [1, 73.5, 147], {'0', '\pi', '2\pi'}, [], {}, 24, [], 0);

%% ---------------------------------
% neuronal firing across phases in vmpfc
% ----------------------------------

curr = figure(7);

confInt = [0, 1];
subplot(6, 6, 2:3)
avgedFR = squeeze(nanmean(nanmean(nanmean(unitPhaseLocked, 2), 1), 4))';
overallAvg = nanmean(avgedFR, 2);

figInfo = {'', 'Phase', 'Firing Rate Change [Hz]', ... 
                [], [], [1, floor(length(auxStruct.radBounds)/2), length(auxStruct.radBounds)], {'0', '\pi', '2\pi'}, [], {}, 24, [], []};

plotBar((avgedFR - nanmean(avgedFR, 2))*1000, figInfo, auxStruct.blueVec, auxStruct.blueVec, curr, confInt, 0);

% check significance with an anova
avgedFR_pop = squeeze(nanmean(nanmean(nanmean(unitPhaseLocked, 2), 1), 4))';

%anovaRes = anova1(avgedFR_pop - nanmean(avgedFR_pop, 2));

%Source      SS       df        MS         F       Prob>F  
%Columns   0.00003      9   3.35181e-06   9.86   8.0192e-15

% repeated measures anove

avg_ph_fr = (avgedFR - nanmean(avgedFR, 2))*1000;
id        = 1:size(avg_ph_fr, 1);
tmp       = table(id', avg_ph_fr(:, 1), avg_ph_fr(:, 2), avg_ph_fr(:, 3), avg_ph_fr(:, 4), avg_ph_fr(:, 5), avg_ph_fr(:, 6), avg_ph_fr(:, 7), avg_ph_fr(:, 8), avg_ph_fr(:, 9), avg_ph_fr(:, 10));
w         = table(categorical([1 2 3 4 5 6 7 8 9 10].'), 'VariableNames', {'phase'}); % within-desin
rm = fitrm(tmp, 'Var2-Var11 ~ 1', 'WithinDesign', w);
ranova(rm, 'withinmodel', 'phase')

%     (Intercept)          8.9055e-30       1    8.9055e-30    9.336e-14             1             1             1            1
%     Error                1.5167e-14     159    9.5389e-17                                                                    
%     (Intercept):phase        30.166       9        3.3518       8.8742    4.2028e-13    2.4897e-05    2.1539e-05    0.0033461
%     Error(phase)             540.49    1431        0.3777                                                                    

%% ---------------------------------
% sample cells showing this modulation, esp over time 
% ----------------------------------
curr = figure(6)

sample_cells = squeeze(nanmean(nanmean(unitPhaseLocked, 1), 2));
i = 27; % sample cell with clear increase before choice + in those phases
subplot(6, 6, 4:5)
imagesc(imgaussfilt(sample_cells(:, :, i)*1000, 1.0)); 
hcb = colorbar;
hcb.Title
hcb.Title.String = "Firing Rate [Hz]";
colormap viridis
figElements(curr, 'Cell Example', 'Time [msec]', 'Phase',... 
                [1 80], [], [30 50 80], {'-500', '-300', 'Choice ON'}, [1, floor(length(auxStruct.radBounds)/2), length(auxStruct.radBounds)-1], {'0', '\pi', '2\pi'}, 24, [], []);  

subplot(6,6,6)
sample_cells_trials = squeeze(nanmean(unitPhaseLocked(:, :, 5, :, 27), 1));
% remove nans 
sample_cells_trials = removenans(sample_cells_trials, 1);
a=plotmse(sample_cells_trials(:, 1:85) * 1000, [230/255, 159/255, 1/255], [0, 1])
sample_cells_trials = squeeze(nanmean(unitPhaseLocked(:, :, 10, :, 27), 1));
sample_cells_trials = removenans(sample_cells_trials, 1);
b=plotmse(sample_cells_trials(:, 1:85) * 1000, [86/255, 180/255, 233/255], [0, 1])
figElements(curr, '', 'Time [msec]', 'Firing Rate [Hz]', [1 80], [], [30 50 80], {'-500', '-300', 'Choice ON'}, [], {}, 24, [], []);  
%legend([a, b], {'Preferred Phase', 'Non-Preferred Phase'}, 'location', 'best')

%% ---------------------------------
% averaged theta oscillation across vmpfc channels
% ----------------------------------

curr = figure(6)
subplot(6, 1, 2)
avgTheta_pl = avgTheta(~isnan(avgTheta(:, 1)), :); % these are vmpfc chans
rectangle('Position', [501, -1.5*10^-3, 300, 3*10^-3], ...
       'FaceColor', [1/255, 114/255, 178/255, 0.2], ...
       'EdgeColor', 'none'); hold on
   
plotmse(avgTheta_pl, auxStruct.grayVec, [0, 1])
figElements(curr, '', 'Time [msec]', 'Amplitude [mV]',...
    [100 1200], [-1.5*10^-3, 1.5*10^-3], [300 500 800], {'-500', '-300', 'Choice onset'}, [], {}, 24, [], []);

%% ---------------------------------
% hex signal across the four peaks 
% ----------------------------------

curr = figure(55);
subPlotIndx = 9:12;
confInt = [0, 1];

for nPos = 1:4
subplot(6, 4, subPlotIndx(nPos))

cells = find(auxStruct.brain_region_cells == 5);
[a,b,c,d]=ttest(stored_betas_roi(cells, :, nPos))

if nPos == 1
    figInfo = {'', 'Symmetries', 'Betas [a.u.]', [], [-0.025, 0.04], [1 2 3 4 5], ...
        {'4', '5', '6', '7', '8'}, [], {}, 24, [], []};
else
    figInfo = {'', 'Symmetries', '', [], [-0.025, 0.04], [1 2 3 4 5], ...
        {'4', '5', '6', '7', '8'}, [], {}, 24, [], 1};
end

[curr] = plotBar(stored_betas_roi(cells, :, nPos), figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, confInt, 0);
title(num2str(nPos));

end

%% ---------------------------------
% example of phase and aligning the phase
% ----------------------------------

curr = figure(5)

offsets = [];
for tr = 1:280
    [a,b] = findpeaks(squeeze(overallPhases(tr, 500:800))', 'MinPeakProminence', 1.5);
    if isempty(b)
        ...
    else
        offsets(tr) = b(end);
    end
end
    
subplot(12, 4,  [25, 29])
imagesc(squeeze(overallPhases(1:201, 500:800))); axis square
colormap viridis
figElements(curr, '', 'Time', 'Trials',... 
                [], [], [], {}, [], {}, 24, 1, 1);  
trIndx = 1;
aligned_data = [];
h=colorbar
h.Label.String = 'Phase [rad]';

for tr = offsets
    
	aligned_data(trIndx, 1:300) = squeeze(overallPhases(trIndx, (tr + 500 - 149) : (tr + 500 + 150)));
    trIndx = trIndx + 1;
end

subplot(12, 4, [26, 30])
imagesc(aligned_data); axis square
colormap viridis
figElements(curr, '', 'Time', '',... 
                [], [], [], {}, [], {}, 24, 1, 1);      
h=colorbar
h.Label.String = 'Phase [rad]';

%% ---------------------------------
% computing + plotting significance of hex signal locked in a phase-aligned way
% ----------------------------------

vmpfc_cells = find(auxStruct.brain_region_cells == 5);
% compute the stats for the empirical grid phase betas
tmp = permute(squeeze(nanmean(gridPhaseBetas(:, 1, :, :, vmpfc_cells), 1)), [3, 1, 2]);
[a,b,c,d]= ttest(tmp);
realStats = d.tstat(:, :, 3);

% compute the states for the null grid phase betas 
permStats = permutationStats(3, :, :);

ou = compute_cluster_stats(realStats(:, auxStruct.neuron_theta_aligned_roi), permStats(:, auxStruct.neuron_theta_aligned_roi, :), 97.5, 3)

%           clustLenUB: 28
%           clustLenLB: 28
%          clustMassUB: 74.4217
%          clustMassLB: 76.4132
%               trueUB: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 … ]
%               trueLB: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 … ]
%            trueLenUB: 32
%            trueLenLB: []
%           trueMassUB: 82.0907
%           trueMassLB: []
%           trueIndxUB: {[32×1 double]}
%           trueIndxLB: {}
%             threshUB: [1.9396 1.9360 1.9638 1.8635 1.8876 1.8547 1.8995 1.9449 1.9778 1.9918 1.9288 1.9456 1.9911 1.9766 1.9501 1.9313 1.9216 1.8498 1.8467 1.8895 … ]
%             threshLB: [2.0707 2.1169 2.0930 2.0665 2.0482 2.0865 2.1476 2.1204 2.1381 2.1043 2.1279 2.1317 2.0887 2.0733 2.0430 1.9975 2.0200 2.0100 2.0151 2.0257 … ]
%      survived_len_UB: [NaN 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 … ]
%      survived_len_LB: NaN
%     survived_mass_UB: [NaN 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 … ]
%     survived_mass_LB: NaN

curr = figure(76)
% --------------------------------
% plot phase aligned hex fig part
%h2=subplot(6, 2, 10)
% plot(1:270, repmat(1.975, [1, 270]), ':k'); hold on
% plot(1:270, repmat(-1.975, [1, 270]), ':k'); hold on

% plot(ou.threshUB, '--k'); hold on

plot([1, 300], [0, 0], ':k', 'linewidth', .6); hold on

avgPhase_pl = squeeze(circ_mean(permute(averagedPhaseInformation(:, :, :), [2, 1, 3])));
avgPhase_pl = squeeze(avgPhase_pl(~sum(avgPhase_pl, 2)==0, :));
plot(normalize_bound(circ_mean(avgPhase_pl(:, auxStruct.neuron_theta_aligned_roi))', -1.95, 4.5), 'k', 'linewidth', 1.5); hold on

tmp = permute(squeeze(nanmean(gridPhaseBetas(:, 1, :, :, vmpfc_cells), 1)), [3, 1, 2]);
[a,b,c,d]= ttest(tmp);
tstats = squeeze(d.tstat);

tstatsTroughOneTime = tstats;

rel_sig = ou.survived_mass_UB;
all_nan_sig = nan(200, 1);
all_nan_sig(find(rel_sig == 1)) = 4.8;

plot(all_nan_sig, 'linewidth', 4, 'color', 'k');

a=plot(tstatsTroughOneTime(auxStruct.neuron_theta_aligned_roi, [1 2 4 5]), 'color', auxStruct.grayVec, 'linewidth', 1.5); hold on; box off 
b=plot(tstatsTroughOneTime(auxStruct.neuron_theta_aligned_roi, [3]), 'color', auxStruct.greenVec, 'linewidth', 4); hold on % 
% legend('4-fold', '5-fold', '7-fold', '8-fold', '6-fold'); legend off
figElements(curr, 'Grid-like Code', 'Time [msec]', 'Hexadirectional Modulation T-value [a.u.]',... 
                [1 270], [-2.5 5], [], {}, [], {}, 24, [], []);  
% legend([a(1), b], {'Control Symmetries', 'Sixfold Symmetry'}, 'location', 'best')            


%% ---------------------------------
% cells with high betas in this analysis
% ----------------------------------

curr = figure(67)

overlap_step                = 1;
smoothing_kernel            = 0.5; % sigma for kernel of 2d traces

phaseGRID_Strength = squeeze(nanmean(tmp(:, auxStruct.grid_peak, :), 2));
[~, idx] = sort(phaseGRID_Strength(:, 3), 'descend');

bo = vmpfc_cells(idx(7));

fr_c1 = squeeze(nanmean(nanmean(averagedFRInformation(:, :, auxStruct.grid_peak, bo), 1), 3)); % 220:251 % 389:412

% compute angle
RegressorContainer = generateAngleRegressor(AuxiliaryCells(bo));

c1_mag  = discretize(AuxiliaryCells(bo).Right_pay, 5);
c1_prob = discretize(AuxiliaryCells(bo).Right_prob, 5);

c2_mag  = discretize(AuxiliaryCells(bo).Left_pay, 5);
c2_prob = discretize(AuxiliaryCells(bo).Left_prob, 5);

% clip the empty values
fr_c1       = zscore(fr_c1(1:length(c1_mag))');
fr_c1_ori   = fr_c1;

% remove trials with identical values
fr_c1_ori     = fr_c1_ori(RegressorContainer.excludeTrials);
c1_mag_o  = c1_mag(RegressorContainer.excludeTrials);
c2_mag_o  = c2_mag(RegressorContainer.excludeTrials);
c1_prob_o = c1_prob(RegressorContainer.excludeTrials);
c2_prob_o = c2_prob(RegressorContainer.excludeTrials);

out = struct;
% compute rate map
out.rate_map.c1.main = compute_rate_map(fr_c1_ori, c1_mag_o, c1_prob_o, 0, overlap_step, smoothing_kernel, 1);
out.rate_map.c2.main = compute_rate_map(fr_c1_ori, c2_mag_o, c2_prob_o, 0, overlap_step, smoothing_kernel, 1);

subplot(6, 4, 17)
imagesc_text(flip(out.rate_map.c1.main.rate_map_avg, 1), flip(out.rate_map.c1.main.rate_map_num, 1)); box off; axis square
caxis([-0.25, 0.25])
hcb=colorbar;
hcb.Title
hcb.Title.String = "FR [z-scored]";
colormap viridis
figElements(curr, 'u_{405} Left', 'Prob', 'Mag', [], [], [1, 5], {'Min', 'Max'}, [1, 5], {'Max', 'Min'}, 24, [], []);

subplot(6, 4, 21)
imagesc_text(flip(out.rate_map.c2.main.rate_map_avg, 1), flip(out.rate_map.c2.main.rate_map_num, 1)); box off; axis square
caxis([-0.25, 0.25])
hcb=colorbar;
hcb.Title
hcb.Title.String = "FR [z-scored]";
colormap viridis
figElements(curr, 'u_{405} Right', 'Prob', 'Mag', [], [], [1, 5], {'Min', 'Max'}, [1, 5], {'Max', 'Min'}, 24, [], []);

bo = vmpfc_cells(idx(5));

fr_c1 = squeeze(nanmean(nanmean(averagedFRInformation(:, :, auxStruct.grid_peak, bo), 1), 3)); % 220:251 % 389:412

% compute angle
RegressorContainer = generateAngleRegressor(AuxiliaryCells(bo));

c1_mag  = discretize(AuxiliaryCells(bo).Right_pay, 5);
c1_prob = discretize(AuxiliaryCells(bo).Right_prob, 5);

c2_mag  = discretize(AuxiliaryCells(bo).Left_pay, 5);
c2_prob = discretize(AuxiliaryCells(bo).Left_prob, 5);

% clip the empty values
fr_c1       = zscore(fr_c1(1:length(c1_mag))');
fr_c1_ori   = fr_c1;

% remove trials with identical values
fr_c1_ori = fr_c1_ori(RegressorContainer.excludeTrials);
c1_mag_o  = c1_mag(RegressorContainer.excludeTrials);
c2_mag_o  = c2_mag(RegressorContainer.excludeTrials);
c1_prob_o = c1_prob(RegressorContainer.excludeTrials);
c2_prob_o = c2_prob(RegressorContainer.excludeTrials);

out = struct;
% compute rate map
out.rate_map.c1.main = compute_rate_map(fr_c1_ori, c1_mag_o, c1_prob_o, 0, overlap_step, smoothing_kernel, 1);
out.rate_map.c2.main = compute_rate_map(fr_c1_ori, c2_mag_o, c2_prob_o, 0, overlap_step, smoothing_kernel, 1);

subplot(6, 4, 18)
imagesc_text(flip(out.rate_map.c1.main.rate_map_avg, 1), flip(out.rate_map.c1.main.rate_map_num, 1)); box off; axis square
 caxis([-0.25, 0.25])
hcb=colorbar;
hcb.Title
hcb.Title.String = "FR [z-scored]";
colormap viridis
figElements(curr, 'u_{643} Left', 'Prob', 'Mag', [], [], [1, 5], {'Min', 'Max'}, [1, 5], {'Max', 'Min'}, 24, [], []);

subplot(6, 4, 22)
imagesc_text(flip(out.rate_map.c2.main.rate_map_avg, 1), flip(out.rate_map.c2.main.rate_map_num, 1)); box off; axis square
 caxis([-0.25, 0.25])
hcb=colorbar;
hcb.Title
hcb.Title.String = "FR [z-scored]";
colormap viridis
figElements(curr, 'u_{643} Right', 'Prob', 'Mag', [], [], [1, 5], {'Min', 'Max'}, [1, 5], {'Max', 'Min'}, 24, [], []);


%% ------------------------
% cpd estimates in a phase-time dependent way
% -------------------------

real_mat = squeeze(nanmean(cpd_inf_perm(:, :, 1)));
null_mat = squeeze(nanmean(cpd_inf_perm(:, :, 2:end)));

out_clu = compute_cluster_stats(real_mat(:, auxStruct.neuron_theta_aligned_roi)', null_mat(auxStruct.neuron_theta_aligned_roi, :), 97.5, 2)

%           clustLenUB: 14.5000
%           clustLenLB: 15
%          clustMassUB: 0.0722
%          clustMassLB: 0
%               trueUB: [301×1 logical]
%               trueLB: [301×1 logical]
%            trueLenUB: [2 12 2 8 3 2 1 29 34 15 1 1 1 2]
%            trueLenLB: [2 14 8]
%           trueMassUB: [0.0099 0.0610 0.0098 0.0394 0.0146 0.0097 0.0049 0.1447 0.1693 0.0751 0.0049 0.0048 0.0049 0.0099]
%           trueMassLB: [-0.0085 -0.0583 -0.0334]
%           trueIndxUB: {1×14 cell}
%           trueIndxLB: {[2×1 double]  [14×1 double]  [8×1 double]}
%             threshUB: [301×1 double]
%             threshLB: [301×1 double]
%      survived_len_UB: [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN … ]
%      survived_len_LB: [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN … ]
%     survived_mass_UB: [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN … ]
%     survived_mass_LB: [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN … ]

length(149:177)
length(213:246)
length(253:267)
find(~isnan(out_clu.survived_len_UB))

cell_cpd = cpd_inf_perm(:, :, 1);

curr = figure(76)
plot(normalize_bound(circ_mean(avgPhase_pl(:, auxStruct.neuron_theta_aligned_roi))', 0.004, 0.0056) * 100, 'k', 'linewidth', 1.5); hold on
plot_sigline_perm(cell_cpd(:, auxStruct.neuron_theta_aligned_roi) * 100, out_clu, 0.55, 1, [0, 1]); box off
%plot(1:301, repmat(prctile(out_clu.threshUB, 97.5)*100, [301, 1]), '--k')
figElements(curr, 'Chosen Value Code', 'Time [msec]', 'CPD [%]',... 
                [1 270], [0.385, .6], [], {}, [], {}, 24, [], []);  

end
