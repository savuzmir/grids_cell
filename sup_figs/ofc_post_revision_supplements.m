function out = run_ofc_post_revision_supplements()
clear; clf;

%%
% cd(go into data folder) 
rel_files = {'unitPhaseLocked_ofc', ...
             'gridPhaseBetas_ofc', 'averagedFRInformation_ofc', 'averagedPhaseInformation_ofc', 'ofc_grid_betas_supp', ...
             'avgTheta_ofc', 'summarised_ofc_cpd_supp', 'auxStruct', 'AuxiliaryCells', 'allUnits'};

for fl = 1:length(rel_files)
    load(rel_files{fl});
end


%% compare amplitude values in ofc to see whether there's a clear peak / trough

ofc = unitPhaseLocked;

confInt = [0, 1];

curr=figure(7)
subplot(2, 2, 1)
avgedFR = squeeze(nanmean(nanmean(nanmean(ofc(:, :, :, :, :), 2), 1), 4));

figInfo = {'', 'Phase', 'Firing Rate Change [Hz]', ... 
                [], [], [1, floor(length(auxStruct.radBounds)/2), length(auxStruct.radBounds)], {'0', '\pi', '2\pi'}, [], {}, 24, [], []};

plotBar((avgedFR - nanmean(avgedFR, 1))'*1000, figInfo, auxStruct.blueVec, auxStruct.blueVec, curr, confInt, 0);

% ofc
avgedFR = squeeze(nanmean(nanmean(nanmean(ofc(:, :, :, :, :), 2), 1), 4));
avg_ph_fr = (avgedFR - nanmean(avgedFR, 1))'*1000;
id        = 1:size(avg_ph_fr, 1);
tmp       = table(id', avg_ph_fr(:, 1), avg_ph_fr(:, 2), avg_ph_fr(:, 3), avg_ph_fr(:, 4), avg_ph_fr(:, 5), avg_ph_fr(:, 6), avg_ph_fr(:, 7), avg_ph_fr(:, 8), avg_ph_fr(:, 9), avg_ph_fr(:, 10));
w         = table(categorical([1 2 3 4 5 6 7 8 9 10].'), 'VariableNames', {'phase'}); % within-desin
rm = fitrm(tmp, 'Var2-Var11 ~ 1', 'WithinDesign', w);
ranova(rm, 'withinmodel', 'phase')


%% get OFC lfp signal 

curr = figure(65)
% plot ofc oscillation
subplot(3, 3, 2:3)
plot(repmat(81, 2, 1), [-.0008, -0.00075], 'linewidth', 5, 'color', 'k'); hold on
ofc_col = auxStruct.areaCol(3, :);
plotmse(q_smooth(avgTheta, 100, 10), ofc_col, [0, 1]);
rectangle('Position', [51, -.0008, 30, .002], ...
'FaceColor', [1/255, 114/255, 178/255, 0.2], ...
'EdgeColor', 'none'); hold on
figElements(curr, '', 'Time [msec]', 'Amplitude [mV]', [1 121], [-.0008, .0008], [41 81 121], ...
    {'-400', 'Choice ON', '400'}, [], {}, 24, [], []); box off

% plot ofc bar plot 
subplot(3, 3, 1)
avgedFR = squeeze(nanmean(nanmean(nanmean(ofc(:, :, :, :, :), 2), 1), 4));

figInfo = {'', 'Phase', 'Firing Rate Change [Hz]', ... 
                [], [], [1, floor(length(auxStruct.radBounds)/2), length(auxStruct.radBounds)], {'0', '\pi', '2\pi'}, [], {}, 24, [], []};

plotBar((avgedFR - nanmean(avgedFR, 1))'*1000, figInfo, auxStruct.blueVec, auxStruct.blueVec, curr, confInt, 0);

% ofc
avgedFR = squeeze(nanmean(nanmean(nanmean(ofc(:, :, :, :, :), 2), 1), 4));
avg_ph_fr = (avgedFR - nanmean(avgedFR, 1))'*1000;
id        = 1:size(avg_ph_fr, 1);
tmp       = table(id', avg_ph_fr(:, 1), avg_ph_fr(:, 2), avg_ph_fr(:, 3), avg_ph_fr(:, 4), avg_ph_fr(:, 5), avg_ph_fr(:, 6), avg_ph_fr(:, 7), avg_ph_fr(:, 8), avg_ph_fr(:, 9), avg_ph_fr(:, 10));
w         = table(categorical([1 2 3 4 5 6 7 8 9 10].'), 'VariableNames', {'phase'}); % within-desin
rm = fitrm(tmp, 'Var2-Var11 ~ 1', 'WithinDesign', w);
ranova(rm, 'withinmodel', 'phase')

%                            SumSq        DF       MeanSq          F           pValue       pValueGG      pValueHF     pValueLB 
%                          __________    ____    __________    __________    __________    __________    __________    _________
% 
%     (Intercept)            1.54e-29       1      1.54e-29    2.9964e-13             1             1             1            1
%     Error                9.9708e-15     194    5.1396e-17                                                                     
%     (Intercept):phase        26.299       9        2.9221        7.8282    2.1772e-11    1.0259e-05    8.5591e-06    0.0056624
%     Error(phase)             651.73    1746       0.37327        
    
figInfo_bl = {'', 'Symmetries', 'Betas [a.u.]', [], [-0.025, 0.04], [1 2 3 4 5], ...
{'4', '5', '6', '7', '8'}, [], {}, 24, [], []};

ofc_col = auxStruct.areaCol(3, :);

out_beta_info = load('ofc_grid_betas_supp.mat');

pl_idx = [5;6;7;8];
for i = 1:4
subplot(3, 4, pl_idx(i, :))
cells = find(auxStruct.brain_region_cells == 3);
[curr] = plotBar(out_beta_info.data_file(cells, :, i), figInfo_bl, ofc_col, ofc_col, curr, confInt, 0); axis square
reshape(print_ttest(out_beta_info.data_file(cells, :, i)), [5, 3])

end

%     0.1655   -1.3919  194.0000
%     0.9809   -0.0240  194.0000
%     0.0435    2.0326  194.0000
%     0.2716    1.1025  194.0000
%     0.1669    1.3875  194.0000


% 
%     0.1642   -1.3963  194.0000
%     0.3722    0.8944  194.0000
%     0.0391    2.0771  194.0000
%     0.0816    1.7506  194.0000
%     0.3133    1.0110  194.0000


%
%     0.0796   -1.7620  194.0000
%     0.8829   -0.1475  194.0000
%     0.1007    1.6493  194.0000
%     0.1849    1.3304  194.0000
%     0.6137    0.5056  194.0000


%     0.0235   -2.2834  194.0000
%     0.1723   -1.3699  194.0000
%     0.2734    1.0985  194.0000
%     0.4223    0.8041  194.0000
%     0.8985   -0.1277  194.0000

        
avgPhase_pl = squeeze(circ_mean(permute(averagedPhaseInformation(:, :, :), [2, 1, 3])));

ofc_ch = find(auxStruct.brain_region_channels == 3);

out_beta_info = load('summarised_ofc_cpd_supp.mat');

curr = figure(7)
% plot cpd value 
subplot(3, 3, 7:8); 
plot(normalize_bound(circ_mean(avgPhase_pl(ofc_ch, auxStruct.neuron_theta_aligned_roi))', 0.45, 0.65), 'k', 'linewidth', 1.5); hold on
real_mat = nanmean(out_beta_info.data_file(:, :, 1), 1);
null_mat = nanmean(out_beta_info.data_file(:, :, 2:end), 1);
out_clu = compute_cluster_stats(real_mat(:, auxStruct.neuron_theta_aligned_roi)', squeeze(null_mat(:, auxStruct.neuron_theta_aligned_roi, 2:end)), 97.5, 2);
%plot(out_clu.threshUB * 100, '--k'); hold on
plot_sigline_perm(out_beta_info.data_file(:, auxStruct.neuron_theta_aligned_roi, 1) * 100, out_clu, 0.63, 1, [0, 1]); box off
figElements(curr, 'Chosen Value Code', 'Time [msec]', 'CPD [%]', [1 300], [0.425, 0.67], [], {}, [], {}, 24, [], []); hold on

% t-test when splitting CPD based on phase pattern 
mn_ph_all = circ_mean(avgPhase_pl(:, auxStruct.neuron_theta_aligned_roi))';
empir_cpd_val = out_beta_info.data_file(:, auxStruct.neuron_theta_aligned_roi, 1);
asc_cpd_val  = squeeze(nanmean(empir_cpd_val(:, find(mn_ph_all>-1&mn_ph_all<1)), 2));
desc_cpd_val = squeeze(nanmean(empir_cpd_val(:, find(mn_ph_all<-1|mn_ph_all>1)), 2));

[a,b,c,d]=ttest(asc_cpd_val, desc_cpd_val)

sig_map =[];
sig_map_p = [];
for i = 1:size(empir_cpd_val, 2)
    for j = 1:size(empir_cpd_val, 2)
        [a,b,c,d] = ttest(empir_cpd_val(:, i, 1), empir_cpd_val(:, j, 1));
        sig_map(i, j) = squeeze(d.tstat);
        sig_map_p(i, j) = b;
    end
end

sig_map(sig_map_p>.05) = 0;

subplot(3, 3, 9); 
imagesc(tril(sig_map, -1)); colorbar; caxis([-3, 3]); colormap viridis; axis square
figElements(curr, 'T-value Map', 'Time [msec]', 'Time [msec]', [], [], [], {}, [], {}, 24, [], []); hold on

%% visualise vmpfc periodic patterns

resid_empir_rate_map = load('vmpfc_supp_rate_maps.mat');
mean_fr = load('vmpfc_supp_rate_maps_mean_fr.mat');

resid_empir_rate_map = resid_empir_rate_map.data_file;
mean_fr              = mean_fr.data_file;

% construct summaries for plotting
rate_maps_summary = nan(5, 5, 9, auxStruct.num_vmpfc_cells);

for i = 1:auxStruct.num_vmpfc_cells

    % left 
    mn = resid_empir_rate_map(:, :, 1, i);
    ub = mn + resid_empir_rate_map(:, :, 2, i);
    lb = mn - resid_empir_rate_map(:, :, 2, i);

    rate_maps_summary(:, :, 1, i) = mn;
    rate_maps_summary(:, :, 2, i) = ub;
    rate_maps_summary(:, :, 3, i) = lb;

    % right
    mn = resid_empir_rate_map(:, :, 4, i);
    ub = mn + resid_empir_rate_map(:, :, 5, i);
    lb = mn - resid_empir_rate_map(:, :, 5);

    rate_maps_summary(:, :, 4, i) = mn;
    rate_maps_summary(:, :, 5, i) = ub;
    rate_maps_summary(:, :, 6, i) = lb;

    % chosen
    mn = resid_empir_rate_map(:, :, 7, i);
    ub = mn + resid_empir_rate_map(:, :, 8, i);
    lb = mn - resid_empir_rate_map(:, :, 8, i);

    rate_maps_summary(:, :, 7, i) = mn;
    rate_maps_summary(:, :, 8, i) = ub;
    rate_maps_summary(:, :, 9, i) = lb;
    i
end

% single cell examples 
cel_vec  = [144, 83, 34, 82, 41, 38, 10, 100, 52];
row_vec  = {[2, 3], [3,4], [2, 3], [1, 2],  [1, 2], [1, 2], [1, 3], [1:5], 2:3};
col_vec  = {[1:5], [1:5], [1:5], [1:5], [1:5], [1:5], [1:5], [2, 5], [1:5]};
indx_vec = {[1:3], [4:6], [1:3], [4:6], [4:6], [1:3], [1:3], [1:3], 1:3};

label = {'Left', 'Right', 'Left', 'Right', 'Right', 'Left', 'Left', 'Left', 'Left'};

curr = figure(6)

pinkVec_a = auxStruct.greenVec + 0.15; % pinkVec - 0.25;
pinkVec_b = auxStruct.greenVec - 0.25; % pinkVec + 0.04;

xvec = repmat((1:5)', [1, 1]);

for i = 1:9
    subplot(3, 3, i)
    plot(1:5, repmat(0, [1, 5]), ':k'); hold on
    cell_inf = ['u', num2str(cel_vec(i)), ' (', num2str(round(mean_fr(cel_vec(i))*1000, 1)), ' Hz)'];
    if i < 8 | i == 9
    a=plotmseSumm(xvec, squeeze(rate_maps_summary(row_vec{i}(1), col_vec{i}, indx_vec{i}, cel_vec(i))), pinkVec_a); box off
    b=plotmseSumm(xvec, squeeze(rate_maps_summary(row_vec{i}(2), col_vec{i}, indx_vec{i}, cel_vec(i))), pinkVec_b); box off
    else
    a=plotmseSumm(xvec, squeeze(rate_maps_summary(row_vec{i}, col_vec{i}(1), indx_vec{i}, cel_vec(i))), pinkVec_a); box off
    b=plotmseSumm(xvec, squeeze(rate_maps_summary(row_vec{i}, col_vec{i}(2), indx_vec{i}, cel_vec(i))), pinkVec_b); box off
    end
    
    if i < 8
        ve_a = ['Mag: ' num2str(row_vec{i}(1))];
        ve_b = ['Mag: ', num2str(row_vec{i}(2))];
    elseif i == 8
        ve_a = ['Prob: ' num2str(col_vec{i}(1))];
        ve_b = ['Prob: ', num2str(col_vec{i}(2))];
    end
    legend([a, b], {ve_a, ve_b}, 'location', 'best'); legend box off 
    
    if i == 1 | i == 4 | i == 7
    figElements(curr, cell_inf, ['Probability_{', label{i}, '}'], 'z-score [a.u.]', [], [-1, 1.5], [1, 5], {'Min', 'Max'}, [], {}, 24, [], [])
    elseif (i > 1 & i < 4) | (i > 4 & i < 8) | (i == 9) 
    figElements(curr, cell_inf, ['Probability_{', label{i}, '}'], '', [], [-1, 1.5], [1, 5], {'Min', 'Max'}, [], {}, 24, [], 1)
    elseif i==8
    figElements(curr, cell_inf, ['Magnitude_{', label{i}, '}'], '', [], [-1, 1.5], [1, 5], {'Min', 'Max'}, [], {}, 24, [], 1)
    end
end
