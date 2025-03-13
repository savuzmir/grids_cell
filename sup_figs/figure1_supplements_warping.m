function out = figure1_supplements_warping(auxStruct)
clear; clf;

%% load data structs
rel_files = {'auxStruct', 'auxLFP', 'AuxiliaryCells', 'chanInfo_regOut', 'subjective_objective_grid', 'allUnits', 'subjective_value_results'};
for fl = 1:length(rel_files)
    load(rel_files{fl});
end

%% grid warping - first fit prospect thoery model
fittedBehaviour = struct;
fittedBehaviour.m1 = [];
fittedBehaviour.m2 = [];
fittedBehaviour.m3 = [];

% selected parameters
% alpha = parameters(1); | x scale
% gamma = parameters(2); | y scale
% theta = parameters(3); | choice temperature

sessionTrials = [];
absIndex = 1;

for an = 1:2

    if an == 1 % need to index correct info
        nSessions = auxStruct.sessionsF;
        trialsSliced = auxStruct.unitSessions(auxStruct.unitSessions(:, 6) == 0, :);
    else
        nSessions = auxStruct.sessionsM;
        trialsSliced = auxStruct.unitSessions(auxStruct.unitSessions(:, 6) == 1, :);
    end

    sessIndx = extract_sess(nSessions, trialsSliced);

    for sess = 1:size(nSessions, 1)
        sessionTrials(absIndex) = size(AuxiliaryCells(sessIndx(sess)).ChML, 1);

        %       mod 1
        startingParam = [auxStruct.alpha, auxStruct.gamma, auxStruct.theta];
        [param, energy] = fit_ProspectTheory1(AuxiliaryCells, sessIndx, sess, startingParam);
        fittedBehaviour.m1(absIndex, :) = [param, energy];

        %       mod 2
        startingParam = [auxStruct.alpha, auxStruct.gamma auxStruct.theta, auxStruct.delta];
        [param, energy] = fit_ProspectTheory2(AuxiliaryCells, sessIndx, sess, startingParam);
        fittedBehaviour.m2(absIndex, :) = [param, energy];

        %       mod 3
        startingParam = [auxStruct.alpha, auxStruct.gamma, auxStruct.theta, auxStruct.delta, auxStruct.zeta];
        [param, energy] = fit_ProspectTheory3(AuxiliaryCells, sessIndx, sess, startingParam);
        fittedBehaviour.m3(absIndex, :) = [param, energy];

        absIndex = absIndex + 1;
    end
end

%% compute bic

logLik = [];
modelBIC = [];

for m = 1:3
    f = {['m', num2str(m)]};
    nParam(m) = size(fittedBehaviour.(f{1}), 2);
    logL = fittedBehaviour.(f{1}); logL = logL(:, end)';
    logLik(m, :) = logL;

    modelBIC(m, :) = log(sessionTrials)*nParam(m) - 2*logL;
end

mean(modelBIC, 2)

% difference
diff(sum(modelBIC, 2))


%% visualise behavioural warping

mags  = normalize_bound(interp1(unique(auxLFP(1).ch.Left_pay)', 1:100, 'spline'), 0.15, 1.5);
probs = normalize_bound(interp1(unique(auxLFP(1).ch.Left_prob)', 1:100, 'spline'), 0.1, 0.9);

spaces = [];
spaces_diff = [];
mag_diff  = [];
prob_diff = [];

sub_mags_tot  = [];
sub_probs_tot = [];

for ss = 1:size(fittedBehaviour.m1, 1)
    out = fittedBehaviour.m1(ss, 1:2);
    tmp_orig = mags .* flip(probs');

    sub_mags  = mags.^out(1);
    sub_probs = exp(-(-log(probs)).^out(2));

    tmp = sub_mags .* flip(sub_probs');
    spaces(:,:,ss)=tmp;
    spaces_diff(:, :, ss) = tmp_orig - tmp;

    mag_diff(:, ss)  = mags - sub_mags;
    prob_diff(:, ss) = probs - sub_probs;

    sub_mags_tot(:, ss)  = sub_mags;
    sub_probs_tot(:, ss) = sub_probs;
end

% plot
curr = figure(7)
surf(nanmean(sub_mags_tot,2), flip(nanmean(sub_probs_tot, 2)), nanmean(spaces_diff, 3), 'EdgeColor', 'none'); hold on
surf(nanmean(sub_mags_tot,2), flip(nanmean(sub_probs_tot, 2)), zeros(100, 100), 'facecolor', [0.1, 0.1, 0.1], 'FaceAlpha', .05, 'EdgeColor', 'k')
figElements(curr, 'EV Difference', 'Magnitude', 'Probability', [0.2, 1.5], [0.1, 1], [], {}, [], {}, 24, [], []);
zlabel('Theoretical-Empirical EV'); colorbar
colormap(auxStruct.curr_palette)

%% fit QF ob + sub and compare whether scaling of behaviour-optimised outperforms original

vals=load('C:\Users\sebas\My Drive\Work\KennerleyLab\Code\Infogathering\cell_paper_revision\master\fin_data\subjective_objective_grid.mat');

curr = figure(7)
originalBeta = squeeze(mean(vals.data_file(1, :, :, :, :, 1), 3));
fittedBeta   = squeeze(mean(vals.data_file(1, :, :, :, :, 2), 3));

original = originalBeta;
fitted   = fittedBeta;

rel_rg = [5];
[unique_fitted_hex, unique_orig_hex, unique_hex_diff, rel_sess] = compute_session_level_stats(auxStruct.animal, auxStruct.session, fitted, original, auxStruct.brain_region_channels, rel_rg, auxStruct.nTR);
hex_diff_values = unique_fitted_hex;

blueCol = [66, 134, 244]./255;
redCol  = [0.6, 0.3, 0.4];
symm    = 3; % corresponds to sixfold

% plot hex effect
subplot(1, 5, 1:2)
a=plotmse(squeeze(unique_orig_hex(symm, auxStruct.rel_sr_tp, find(rel_sess>auxStruct.trialLim)))', blueCol, [0, 1]); hold on
b=plotmse(squeeze(unique_fitted_hex(symm, auxStruct.rel_sr_tp, find(rel_sess>auxStruct.trialLim)))', redCol, [0, 1])
plot([1, 200], [0, 0], ':k', 'linewidth', 0.6);
legend([a, b], {'Objective', 'Subjective'})
figElements(curr, '', 'Time [msec]', 'Betas [a.u.]', [1, 121], [-.1, .15], [21, 61, 101], ...
    {'-400', 'Choice ON', '400'}, [], {}, 24, [], []); box off

% stats
inps = [];
inps(:, :, 1) = squeeze(hex_diff_values(symm, auxStruct.rel_sr_tp, :))';
inps(:, :, 2) = zeros(54, length(auxStruct.rel_sr_tp)); % cluster against 0
nan_rw = find(isnan(inps(:, 1, 1)));
inps(nan_rw, :, 2) = nan;
inps(find(rel_sess<auxStruct.trialLim), :, :) = nan;
rng(07062023)
ou = run_cluster_based_perm_test(inps, 1, auxStruct.nPerm, 1, 95, 2);
ou.cluster_stats
plot(ou.cluster_stats.survived_mass_UB - 0.9, 'color', blueCol, 'linewidth', 2)
legend([a, b], {'Objective Sixfold Symmetry', 'Subjective Sixfold Symmetry'}, 'box', 'off')

% plot indiv plots for these regions in original roi
% ====== in vmpfc
confInt = [0, 0];
pl_idx = [3, 4, 5];
rg_idx = [5, 1, 3];
nam_idx = {'VMPFC', 'ACC', 'OFC'};
for rg_l = 1:3
    rg = rg_idx(rg_l);
    pl = pl_idx(rg_l);
    rel_rg = rg;
    [unique_fitted_hex, unique_orig_hex, ~, rel_sess] = compute_session_level_stats(auxStruct.animal, auxStruct.session, fitted, original, auxStruct.brain_region_channels, rel_rg, auxStruct.nTR);
    rel_rws = find(rel_sess > auxStruct.trialLim);

    subplot(1, 5, pl)
    figInfo = {nam_idx{rg_l}, 'Grid-like Code', '', [], [-0.015, 0.1], [1 2], ...
        {'Obj.', 'Subj.'}, [], {}, 24, [], []};
    [curr] = plotBar([squeeze(nanmean(unique_orig_hex(3, auxStruct.inputT, rel_rws), 2)), squeeze(nanmean(unique_fitted_hex(3, auxStruct.inputT, rel_rws), 2))], figInfo, blueCol, redCol, curr, confInt, 0);
    print_ttest(squeeze(nanmean(unique_orig_hex(3, auxStruct.inputT, rel_rws), 2)), 1, squeeze(nanmean(unique_fitted_hex(3, auxStruct.inputT, rel_rws), 2)))
    % 0.0253    2.5036   14.0000
end

%% subjective objective model fitting
% fit models to behaviour

[a,b,c,d]=ttest(modelBIC(1, :), modelBIC(2, :))
[a,b,c,d]=ttest(modelBIC(2, :), modelBIC(3, :))

% subjective much better in behaviour
print_ttest(data_file.beh_pred(:, 2), 1, data_file.beh_pred(:, 1))
% 5.06977121613628e-10          7.59023681382699                        53
[a,b,c,d]=ttest(data_file.beh_pred(:, 1), data_file.beh_pred(:, 3))

orangeVec = [1.0000    0.7569    0.0275];
purpleVec = [0.6235    0.2902    0.5882];

curr = figure(7)
figInfo = {'', 'EV Difference', 'Session Model Fit [beta]', ...
    [0, 3], [], [1, 2], {'Objective', 'Subjective'}, [], {}, 24, [], []};
plotBar(data_file.beh_pred(:, 1:2), figInfo, orangeVec, purpleVec, curr, [0, 0], 1, 1); axis square

curr = figure(9)
% compare LR
concat_ob     = [];
concat_sub    = [];
concat_sub_mag  = [];
concat_sub_prob = [];

for i = 1:54
    [unqVals, tmp]    = unique(round(data_file.ob_estimates_lr{i}, 2));
    unq_sub_estimates = round(data_file.sub_estimates_lr{i}, 2);
    unq_sub_estimates = unq_sub_estimates(tmp);

    concat_ob = [concat_ob; unqVals];
    concat_sub = [concat_sub; unq_sub_estimates];

    % repeat for mag
    unq_sub_estimates = data_file.sub_estimates_mag{i};
    concat_sub_mag = [concat_sub_mag; unq_sub_estimates];

    % repeat for prob
    unq_sub_estimates = data_file.sub_estimates_prob{i};
    concat_sub_prob = [concat_sub_prob; unq_sub_estimates];
end

[~, fo] = sort(concat_ob);

overall_sub_sr = concat_sub(fo);
overall_ob_sr  = concat_ob(fo);
zoo = state_mean(round(abs(overall_ob_sr), 2), round(abs(overall_sub_sr), 2));
scatter(zoo(2:end, 2), zoo(2:end, 1), 20, '', 'markerEdgeColor', 'k', 'linewidth', 2); hold on % we don't plot equal difficulty
for i = 2:size(zoo, 1)
    x_val = [zoo(i, 2) - zoo(i, 4), zoo(i, 2) + zoo(i, 4)];
    y_val = [zoo(i, 1), zoo(i, 1)];
    line(x_val, y_val, 'color', auxStruct.grayVec, 'linewidth', 3);
end
plot([0, .85], [0, .85], '--k')
figElements(curr, 'EV Difference', 'Subjective', 'Objective', [0, 0.86], [0, 0.86], [], {}, [], {}, 24, [], [])

nanmean(fittedBehaviour.m1, 1)

curr = figure(7)
mag_val_ob  = unique(AuxiliaryCells(1).Left_pay);
[zoo] = state_mean(repmat(mag_val_ob, [54, 1]), concat_sub_mag);
subplot(2, 1, 1)
% mag separately
scatter(zoo(:, 2), zoo(:, 1), 50, '', 'markerEdgeColor', 'k', 'linewidth', 2); hold on
for i = 1:size(zoo, 1)
    x_val = [zoo(i, 2) - zoo(i, 4), zoo(i, 2) + zoo(i, 4)];
    y_val = [zoo(i, 1), zoo(i, 1)];
    line(x_val, y_val, 'color', auxStruct.grayVec, 'linewidth', 3);
end
plot([0, 1], [0, 1], '--k')
figElements(curr, 'Reward Magnitude', 'Subjective', 'Objective', [0, 1], [0, 1], [], {}, [], {}, 24, 1, 1)

prob_val_ob  = unique(AuxiliaryCells(1).Left_prob);
[zoo] = state_mean(repmat(prob_val_ob, [54, 1]), concat_sub_prob);

subplot(2, 1, 2)
scatter(zoo(:, 2), zoo(:, 1), 50, '', 'markerEdgeColor', 'k', 'linewidth', 2); hold on
for i = 1:size(zoo, 1)
    x_val = [zoo(i, 2) - zoo(i, 4), zoo(i, 2) + zoo(i, 4)];
    y_val = [zoo(i, 1), zoo(i, 1)];
    line(x_val, y_val, 'color', auxStruct.grayVec, 'linewidth', 3);
end
plot([0, 1], [0, 1], '--k')
figElements(curr, 'Reward Probability', 'Subjective', 'Objective', [0, 1], [0, 1], [], {}, [], {}, 24, 1, 1);

%% subjective / objective value code in dlpfc
stat_info = [];
indiv_reg = nan(300, 4);

for i = 1:4

    % ch val
    if i ~= 4
        a = squeeze(data_file.neural_signal(auxStruct.brain_region_cells == auxStruct.areaCodes(i), 1, 1));
        b = squeeze(data_file.neural_signal(auxStruct.brain_region_cells == auxStruct.areaCodes(i), 1, 2));
    else
        a = squeeze(nanmean(data_file.neural_signal(auxStruct.brain_region_cells == auxStruct.areaCodes(i), 1:3, 1), 2));
        b = squeeze(nanmean(data_file.neural_signal(auxStruct.brain_region_cells == auxStruct.areaCodes(i), 1:3, 2), 2));
    end

    indiv_reg(1:size(a, 1), i) = b - a;

    [a,b,c,d]=ttest(b-a)
    stat_info(i, :) = [d.tstat, d.df, b];
end

%          -0.25443861562117                       197         0.799422115421869
%           3.14862779342283                       155       0.00196823809731854
%           2.02955232876639                       194        0.0437679881921433
%          0.798331153255466                       159         0.425869126674725

curr = figure(7)
b = bar(nanmean(indiv_reg, 1), 'edgeColor', 'none'); hold on; box off
b.FaceColor = 'flat';
for j = 1:length(auxStruct.areaCol)
    b.CData(j, :) = auxStruct.areaCol(j, :);
end
err = nanstd(indiv_reg, [], 1)./sqrt(sum(~isnan(indiv_reg)));
errorbar(1:4, nanmean(indiv_reg, 1), err,  '.', 'CapSize', 0, 'linewidth', 2, 'color', [155/255, 155/255, 155/255])
figElements(curr, 'Chosen Value_{Subj.-Obj.}', 'Brain Regions', 'CPD [%]', [], [], [1 2 3 4], {'ACC', 'DLPFC', 'OFC', 'VMPFC'}, [], {}, 24, [], [])

% post-revision info additions

%% post-revision additions
% num channels and cells

% ------------
% num cells + num channels 
cell_num = [];
chan_num = [];
for i = auxStruct.areaCodes
    cell_num(i) = sum(auxStruct.brain_region_cells == i);
    chan_num(i) = sum(auxStruct.brain_region_channels == i);
end

% 198   156   195     0   160
% 149   104   133     0    97

%---------------------
%reaction time for fig1
sessRT_avg = [];
for i = 1:length(auxLFP)
    sessRT_avg(i) = nanmean(auxStruct.sessRT{i});
end

nanmean(sessRT_avg);


%% samples of
tmp = load('C:\Users\sebas\My Drive\Work\KennerleyLab\Code\Infogathering\cell_paper_revision\master\fin_data\chanInfo_regOut.mat');
chanInfo_regOut = tmp.out_new;  

rate_maps_left   = nan(5, 5, 497);
rate_maps_right  = nan(5, 5, 497);
rate_maps_chosen = nan(5, 5, 497);

for i = 1:length(auxStruct.brain_region_channels)
    left_mag   = discretize(chanInfo_regOut(i).Aux.Left_pay, 5);
    right_mag  = discretize(chanInfo_regOut(i).Aux.Right_pay, 5);
    left_prob  = discretize(chanInfo_regOut(i).Aux.Left_prob, 5);
    right_prob = discretize(chanInfo_regOut(i).Aux.Right_prob, 5);

    chosen_mag  = discretize(chanInfo_regOut(i).Aux.ChML + chanInfo_regOut(i).Aux.ChMR, 5);
    chosen_prob = discretize(chanInfo_regOut(i).Aux.ChPL + chanInfo_regOut(i).Aux.ChPR, 5);

    % realign to keep consistent with regout
    x1 = chanInfo_regOut(i).Aux.Right_pay - chanInfo_regOut(i).Aux.Left_pay;
    y1 = chanInfo_regOut(i).Aux.Right_prob - chanInfo_regOut(i).Aux.Left_prob;

    angles = atan2(y1, x1)*180/pi;
    [~, tmp] = sort(mod(angles, 360));

    curr_val = nanmean(chanInfo_regOut(i).Y(:, auxStruct.PE_roi), 2);

    % compute rate map
    out_left  = compute_rate_map(curr_val, left_mag(tmp), left_prob(tmp), 0, 1, 0, 1);
    out_right = compute_rate_map(curr_val, right_mag(tmp), right_prob(tmp), 0, 1, 0, 1);
    out_chosen = compute_rate_map(curr_val, chosen_mag(tmp), chosen_prob(tmp), 0, 1, 0, 1);

    rate_maps_left(:, :, i)  = out_left.rate_map_avg;
    rate_maps_right(:, :, i) = out_right.rate_map_avg;
    rate_maps_chosen(:, :, i) = out_chosen.rate_map_avg;

    rate_maps_left_s(:, :, i)  = out_left.rate_map_sem;
    rate_maps_right_s(:, :, i) = out_right.rate_map_sem;
    rate_maps_chosen_s(:, :, i) = out_chosen.rate_map_sem;
    i
end

% construct summaries for plotting
rate_maps_summary = [];

for i = auxStruct.relRegs{4}


    mn = rate_maps_left(:, :, i);
    ub = mn + rate_maps_left_s(:, :, i);
    lb = mn - rate_maps_left_s(:, :, i);

    % left
    rate_maps_summary(:, :, 1, i) = mn;
    rate_maps_summary(:, :, 2, i) = ub;
    rate_maps_summary(:, :, 3, i) = lb;

    mn = rate_maps_right(:, :, i);
    ub = mn + rate_maps_right_s(:, :, i);
    lb = mn - rate_maps_right_s(:, :, i);

    % right
    rate_maps_summary(:, :, 4, i) = mn;
    rate_maps_summary(:, :, 5, i) = ub;
    rate_maps_summary(:, :, 6, i) = lb;

end

% lfp examples
% prob (rows) x mag (cols) map

cel_vec  = [ 252,  175,   234,   434,   166,   234,   232,   250,   253,   445,   449,   234,   406];
row_vec  = {[3,5], [1,2], [1,2], [1, 4], [1,3], [3,5], [3,4], [1,2], [1:5], [1:5], [1:5], [1:5], [1:5]}; % mags
col_vec  = {[1:5], [1:5], [1:5], [1:5], [1:5], [1:5], [1:5], [1:5],  [4,5], [1,2], [3,5], [2,3], [2,4]}; % probs
indx_vec = {[1:3], [1:3], [1:3], [4:6], [4:6], [1:3], [4:6], 1:3,   [4:6], [4:6], [4:6], 1:3,   1:3}; % left  / right
label = {'Left', 'Left', 'Left', 'Right', 'Right', 'Left', 'Right', 'Left', 'Right', 'Right', 'Right', 'Left', 'Left'};

corr_indx = [2, 5, 3, 4, 8, 13, 7, 11, 9];
pinkVec_a = auxStruct.greenVec + 0.15; % pinkVec - 0.25;
pinkVec_b = auxStruct.greenVec - 0.25; % pinkVec + 0.04;

curr = figure(6)
xvec = repmat((1:5)', [1, 1]);

pic_indx = 1;

for i = corr_indx
    subplot(3, 3, pic_indx)
    plot(1:5, repmat(0, [1, 5]), ':k'); hold on
    cell_inf = ['ch', num2str(cel_vec(i))];
    if length(col_vec{i}) > 2 %
        a=plotmseSumm(xvec, squeeze(rate_maps_summary(row_vec{i}(1), col_vec{i}, indx_vec{i}, cel_vec(i))), pinkVec_a); box off
        b=plotmseSumm(xvec, squeeze(rate_maps_summary(row_vec{i}(2), col_vec{i}, indx_vec{i}, cel_vec(i))), pinkVec_b); box off
    else
        a=plotmseSumm(xvec, squeeze(rate_maps_summary(row_vec{i}, col_vec{i}(1), indx_vec{i}, cel_vec(i))), pinkVec_a); box off
        b=plotmseSumm(xvec, squeeze(rate_maps_summary(row_vec{i}, col_vec{i}(2), indx_vec{i}, cel_vec(i))), pinkVec_b); box off;
    end

    if length(col_vec{i}) > 2
        ve_a = ['Mag: ' num2str(row_vec{i}(1))];
        ve_b = ['Mag: ', num2str(row_vec{i}(2))];
    else
        ve_a = ['Prob: ' num2str(col_vec{i}(1))];
        ve_b = ['Prob: ', num2str(col_vec{i}(2))];
    end

    legend([a, b], {ve_a, ve_b}, 'location', 'best', 'orientation', 'horizontal'); legend box off

    if pic_indx == 1 | pic_indx == 4 | pic_indx == 7
        y_lab = 'z-score [a.u.]';
        rm_y = [];
    else
        y_lab = '';
        rm_y = 1;
    end

    y_lab = 'z-score [a.u.]';
    rm_y = [];

    if length(col_vec{i}) > 2
        figElements(curr, cell_inf, ['Probability_{', label{i}, '}'], y_lab, [1 5], [-1, 1.5], [1, 5], {'Min', 'Max'}, [], {}, 24, [], rm_y)
    else
        figElements(curr, cell_inf, ['Magnitude_{', label{i}, '}'], y_lab, [1 5], [-1, 1.5], [1, 5], {'Min', 'Max'}, [], {}, 24, [], rm_y)
    end

    pic_indx = pic_indx + 1;
end

