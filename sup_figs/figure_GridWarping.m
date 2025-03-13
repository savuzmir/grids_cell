function figure_GridWarping()

% fit behaviour 

% fit models to behaviour

fittedBehaviour = struct;
fittedBehavour.m1 = [];
fittedBehavour.m2 = [];
fittedBehavour.m3 = [];

% selected parameters
% 
% alpha = parameters(1); | x scale
% gamma = parameters(2); | y scale
% theta = parameters(3); | choice temperature
% delta = parameters(4); | noise parameter
% eta   = parameters(5); | integration coeff
% beta  = parameters(6); | weighing of mag/prob
% zeta  = parameters(7); | side bias

alpha   = 1;
gamma   = 1;
theta   = 1;
absIndex = 1;
sessionTrials = [];

for an = 1:2
    
    if an == 1
        nSessions = auxStruct.sessionsF;
        trialsSliced = auxStruct.unitSessions(auxStruct.unitSessions(:, 6) == 0, :);
    else
        nSessions = auxStruct.sessionsM;
        trialsSliced = auxStruct.unitSessions(auxStruct.unitSessions(:, 6) == 1, :);
    end

    sessIndx = extract_sess(nSessions, trialsSliced);

    for sess = 1:size(nSessions, 1)

        sessionTrials(absIndex) = size(AuxiliaryCells(sessIndx(sess)).ChML, 1);
        
%         mod 1: alpha, gamma, theta
        startingParam = [alpha, gamma, theta]; % just alpha, gamma, theta
        [param, energy] = fit_ProspectTheory1(AuxiliaryCells, sessIndx, sess, startingParam);
        fittedBehaviour.m1(absIndex, :) = [param, energy];
        
%         mod 2: alpha, gamma, theta, delta
        startingParam = [alpha, gamma theta, delta]; % alpha, gamma, theta, delta
        [param, energy] = fit_ProspectTheory2(AuxiliaryCells, sessIndx, sess, startingParam);
        fittedBehaviour.m2(absIndex, :) = [param, energy];
%         
%         mod 3: alpha, gamma, theta, delta, zeta
        startingParam = [alpha, gamma, theta, delta, zeta]; % alpha gamma, theta, delta, zeta
        [param, energy] = fit_ProspectTheory3(AuxiliaryCells, sessIndx, sess, startingParam);
        fittedBehaviour.m3(absIndex, :) = [param, energy];
        
        absIndex = absIndex + 1;
    end
end

% compute bic 

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

beh_pred = [];
neural_signal = nan(724, 1, 4);
ss = 1;

% also test vmPFC within subjective window
load('C:\Users\sebas\OneDrive - University College London\Sebastijan\KennerleyLab\Data\FIL_update\simultaneous_trials\cells_new\averagedFRInformation.mat');

vmpfc_indx = 1;

ob_estimates  = {};
sub_estimates = {};

for unit = 1:size(allUnits, 2)
    if unit > 1

        curr_sess = auxStruct.unitSessions(unit, 2);
        if curr_sess ~= prev_sess
            ss = ss + 1;
            prev_sess = curr_sess;
        end

    else
        curr_sess = auxStruct.unitSessions(unit, 2);
        prev_sess = curr_sess;
    end

    foo = fittedBehaviour.m1(ss, 1:2);
  
    % define objective info 
    left_val  = [AuxiliaryCells(unit).Left_pay .* AuxiliaryCells(unit).Left_prob]; 
    right_val = [AuxiliaryCells(unit).Right_pay .* AuxiliaryCells(unit).Right_prob];

    ob_estimates{ss} = abs(left_val - right_val);

    chDiff   = right_val-left_val;
    choice   = AuxiliaryCells(unit).ChML == 0;

    b = glmfit(chDiff, choice, 'binomial','Link','probit');
    beh_pred(ss, 1) = b(2);

    % subjective
    left_val  = AuxiliaryCells(unit).Left_pay.^foo(1) .* exp(-(-log(AuxiliaryCells(unit).Left_prob)).^foo(2)); 
    right_val = AuxiliaryCells(unit).Right_pay.^foo(1) .* exp(-(-log(AuxiliaryCells(unit).Right_prob)).^foo(2)); 
    
    chDiff  = right_val-left_val;
    b = glmfit(chDiff, choice, 'binomial','Link','probit');
    beh_pred(ss, 2) = b(2);
    
    % frontal area coding, e.g. chosen value and ch val diff 
    chosen_value_ob  = (AuxiliaryCells(unit).ChML + AuxiliaryCells(unit).ChMR) .* (AuxiliaryCells(unit).ChPL + AuxiliaryCells(unit).ChPR);
    chosen_value_sub = ((AuxiliaryCells(unit).ChML + AuxiliaryCells(unit).ChMR).^foo(1)) .* exp(-(-log((AuxiliaryCells(unit).ChPL + AuxiliaryCells(unit).ChPR))).^foo(2));
    
    ob_estimates{ss} = chosen_value_ob; 
    sub_estimates{ss} = chosen_value_sub;

    unchosen_value_ob  = (AuxiliaryCells(unit).UnML + AuxiliaryCells(unit).UnMR) .* (AuxiliaryCells(unit).UnPL + AuxiliaryCells(unit).UnPR);
    unchosen_value_sub = ((AuxiliaryCells(unit).UnML + AuxiliaryCells(unit).UnMR).^foo(1)) .* exp(-(-log((AuxiliaryCells(unit).UnPL + AuxiliaryCells(unit).UnPR))).^foo(2));
        
    rel_y_dat = allUnits(unit).trials;
    full_dm_ob   = [ones(size(rel_y_dat, 1), 1), chosen_value_ob];
    full_dm_sub  = [ones(size(rel_y_dat, 1), 1), chosen_value_sub];

    full_dm_ob_diff  = [ones(size(rel_y_dat, 1), 1), chosen_value_ob - unchosen_value_ob];
    full_dm_sub_diff = [ones(size(rel_y_dat, 1), 1), chosen_value_sub - unchosen_value_sub];

    if brain_regionCells(unit) == 5
        rel_y_dat = squeeze(averagedFRInformation(:, :, :, unit));
        vmpfc_indx = vmpfc_indx + 1;

        tmp_a = [];
        tmp_b = [];

        for fr = 1:14

                curr_rel_y_dat = squeeze(rel_y_dat(fr, 1:size(full_dm_ob, 1), :));

                [a]     = cpd(nanmean(curr_rel_y_dat(:, 245), 2), full_dm_ob); % test on ROIs
                [b]     = cpd(nanmean(curr_rel_y_dat(:, 245), 2), full_dm_sub);
                tmp_a(fr, 1, :) = a(2, :);
                tmp_b(fr, 1, :) = b(2, :);                 
                [a]     = cpd(nanmean(curr_rel_y_dat(:, 343), 2), full_dm_ob); % test on ROIs
                [b]     = cpd(nanmean(curr_rel_y_dat(:, 343), 2), full_dm_sub);
                tmp_a(fr, 2, :) = a(2, :);
                tmp_b(fr, 2, :) = b(2, :); 
                [a]     = cpd(nanmean(curr_rel_y_dat(:, 392), 2), full_dm_ob); % test on ROIs
                [b]     = cpd(nanmean(curr_rel_y_dat(:, 392), 2), full_dm_sub);
                tmp_a(fr, 3, :) = a(2, :);
                tmp_b(fr, 3, :) = b(2, :); 

                % difference
                [a]     = cpd(nanmean(curr_rel_y_dat(:, 245), 2), full_dm_ob_diff); % test on ROIs
                [b]     = cpd(nanmean(curr_rel_y_dat(:, 245), 2), full_dm_sub_diff);
                tmp_a(fr, 4, :) = a(2, :);
                tmp_b(fr, 4, :) = b(2, :);                 
                [a]     = cpd(nanmean(curr_rel_y_dat(:, 343), 2), full_dm_ob_diff); % test on ROIs
                [b]     = cpd(nanmean(curr_rel_y_dat(:, 343), 2), full_dm_sub_diff);
                tmp_a(fr, 5, :) = a(2, :);
                tmp_b(fr, 5, :) = b(2, :); 
                [a]     = cpd(nanmean(curr_rel_y_dat(:, 392), 2), full_dm_ob_diff); % test on ROIs
                [b]     = cpd(nanmean(curr_rel_y_dat(:, 392), 2), full_dm_sub_diff);
                tmp_a(fr, 6, :) = a(2, :);
                tmp_b(fr, 6, :) = b(2, :);                 
 
        end

        a = squeeze(nanmean(tmp_a(:, :, :), 1));
        b = squeeze(nanmean(tmp_b(:, :, :), 1));

        a_len = 6;

    else

        [a1]     = cpd(nanmean(rel_y_dat(:, 51:81), 2), full_dm_ob);  a1 = a1(2, :);
        [b1]     = cpd(nanmean(rel_y_dat(:, 51:81), 2), full_dm_sub); b1 = b1(2, :);    

        [a2]     = cpd(nanmean(rel_y_dat(:, 51:81), 2), full_dm_ob_diff);  a2 = a2(2, :);
        [b2]     = cpd(nanmean(rel_y_dat(:, 51:81), 2), full_dm_sub_diff); b2 = b2(2, :);    

        a = [a1;a2];
        b = [b1;b2];

        a_len  = 2;
    end
        neural_signal(unit, 1:a_len, 1) = a'; % objective
        neural_signal(unit, 1:a_len, 2) = b'; % subjective
unit
end



% subjective much better in behaviour
[a,b,c,d]=ttest(beh_pred(:, 2), beh_pred(:, 1))

orangeVec = [1.0000    0.7569    0.0275];
purpleVec = [0.6235    0.2902    0.5882];

curr=figure(6)
subplot(2, 2, 1)
figInfo = {'', 'EV Difference', 'Session Model Fit [GLM]', ... 
                [], [], [1, 2], {'Objective', 'Subjective'}, [], {}, 24, [], []};
plotBar(beh_pred, figInfo, orangeVec, purpleVec, curr, confInt, 1, 1);

[a,b,c,d]=ttest(beh_pred(1:29, 1), beh_pred(1:29, 2))
[a,b,c,d]=ttest(beh_pred(30:54, 1), beh_pred(30:54, 2))

stat_info = [];
indiv_reg = nan(300, 4);

for i = 1:4

    % ch val
    if i ~= 4
        a = squeeze(neural_signal(brain_regionCells == areaCodes(i), 1, 1));
        b = squeeze(neural_signal(brain_regionCells == areaCodes(i), 1, 2));
    else
        a = squeeze(nanmean(neural_signal(brain_regionCells == areaCodes(i), 1:3, 1), 2));
        b = squeeze(nanmean(neural_signal(brain_regionCells == areaCodes(i), 1:3, 2), 2));
    end

    indiv_reg(1:size(a, 1), i) = b - a;

    [a,b,c,d]=ttest(b-a)
    stat_info(i, :) = [d.tstat, b];
end

subplot(2, 2, 2)
b = bar(nanmean(indiv_reg, 1), 'edgeColor', 'none'); hold on; box off 
b.FaceColor = 'flat';
for j = 1:length(areaCol)
    b.CData(j, :) = areaCol(j, :);
end
err = nanstd(indiv_reg, [], 1)./sqrt(sum(~isnan(indiv_reg)));
errorbar(1:4, nanmean(indiv_reg, 1), err,  '.', 'CapSize', 0, 'linewidth', 2, 'color', [155/255, 155/255, 155/255])
figElements(curr, 'Chosen Value', 'Brain Regions', 'CPD [%]', [], [], [1 2 3 4], {'ACC', 'DLPFC', 'OFC', 'VMPFC'}, [], {}, 24, [], [])

stat_info = [];
indiv_reg = nan(300, 4);

for i = 1:4
    % ch val diff
    if i ~= 4
        a = squeeze(neural_signal(brain_regionCells == areaCodes(i), 2, 1));
        b = squeeze(neural_signal(brain_regionCells == areaCodes(i), 2, 2));
    else
        a = squeeze(nanmean(neural_signal(brain_regionCells == areaCodes(i), 4:6, 1), 2));
        b = squeeze(nanmean(neural_signal(brain_regionCells == areaCodes(i), 4:6, 2), 2));
    end
    indiv_reg(1:size(a, 1), i) = b - a;

    [a,b,c,d]=ttest(b-a)
    stat_info(i, :) = [d.tstat, b];
end



% dlpfc ch vs. chdiff fit
a = squeeze(neural_signal(brain_regionCells == 2, 1, 1));
b = squeeze(neural_signal(brain_regionCells == 2, 1, 2));
ch_dlpfc = b-a;

a = squeeze(neural_signal(brain_regionCells == 2, 2, 1));
b = squeeze(neural_signal(brain_regionCells == 2, 2, 2));
ch_diff_dlpfc = b-a;

[a,b,c,d]=ttest(ch_diff_dlpfc, ch_dlpfc);

subplot(2, 2, 3)
b = bar(nanmean(indiv_reg, 1), 'edgeColor', 'none'); hold on; box off 
b.FaceColor = 'flat';
for j = 1:length(areaCol)
    b.CData(j, :) = areaCol(j, :);
end
err = nanstd(indiv_reg, [], 1)./sqrt(sum(~isnan(indiv_reg)));
errorbar(1:4, nanmean(indiv_reg, 1), err,  '.', 'CapSize', 0, 'linewidth', 2, 'color', [155/255, 155/255, 155/255])
figElements(curr, 'Chosen Value Difference', 'Brain Regions', 'CPD [%]', [], [], [1 2 3 4], {'ACC', 'DLPFC', 'OFC', 'VMPFC'}, [], {}, 24, [], [])

% 

overall_sub = [];
overall_ob  = [];
for i = 1:54
    overall_sub = [overall_sub; sub_estimates{i}];
    overall_ob  = [overall_ob; ob_estimates{i}];
end

overall_sub_sr = overall_sub(fo);
overall_ob_sr  = overall_ob(fo);
subplot(2, 2, 4)
zoo = state_mean(overall_ob_sr, overall_sub_sr)
scatter(zoo(:, 2), zoo(:, 1), 20, '', 'markerEdgeColor', 'k', 'linewidth', 2); hold on
for i = 1:size(zoo, 1)
    x_val = [zoo(i, 2) - zoo(i, 4), zoo(i, 2) + zoo(i, 4)];
    y_val = [zoo(i, 1), zoo(i, 1)];
    line(x_val, y_val, 'color', grayVec, 'linewidth', 3);
end
plot([0, 1], [0, 1], '--k')
figElements(curr, 'Chosen Value', 'Subjective', 'Objective', [0, 0.86], [0, 0.86], [], {}, [], {}, 24, [], [])
saveImage(curr, ['objective_subjective'], 'C:\Users\sebas\OneDrive - University College London\Sebastijan\KennerleyLab\Figures\grids_working\paper_figures\revision\', [12000, 12000])



curr = figure(7);


spaces = [];
spaces_diff = [];
mag_diff  = [];
prob_diff = [];

mags = 0.001:0.2:1;
probs = 0.001:0.2:1;

sub_mags_tot  = [];
sub_probs_tot = [];

for ss = 1:54
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

subplot(2, 3, 1)
surf(nanmean(sub_mags_tot,2), flip(nanmean(sub_probs_tot, 2)), nanmean(spaces_diff, 3), 'EdgeColor', 'none'); hold on
surf(nanmean(sub_mags_tot,2), flip(nanmean(sub_probs_tot, 2)), zeros(100, 100), 'facecolor', [0.1, 0.1, 0.1], 'FaceAlpha', .05, 'EdgeColor', 'k')
figElements(curr, 'EV Difference', 'Magnitude States', 'Probability States', [0.2, 1.5], [0.1, 1], [], {}, [], {}, 18, [], []);
zlabel('Theoretical-Empirical EV'); colorbar
colormap(curr_palette)


subplot(2, 5, 6:7) % we then plotted this over time, as in Figure 2B
rel_rg = [1 3 5];
symm = 3;
[unique_fitted_hex, unique_orig_hex, unique_hex_diff, rel_sess] = compute_session_level_stats(animal, session, fitted, original,brain_region, rel_rg, nTR);
rel_rws = find(rel_sess > 200);
plot(1:150, repmat(0, [1, 150]), '--k')
a=plotmse(squeeze(unique_orig_hex(symm, 1:end, rel_rws))', blueCol, [0, 1]); hold on
b=plotmse(squeeze(unique_fitted_hex(symm, 1:end, rel_rws))', redCol, [0, 1])
figElements(curr, 'Veridical and Distorted Hexadirectional Modulation', 'Time [msec]', 'Betas [a.u.]', [20, 150], [], [41, 81, 121], ...
                {'-400', 'Choice ON', '+400'}, [], {}, 14, [], []); box off
axis square 

% there is significantly stronger veridical before and distorted after choice
inps = [];
inps(:, :, 1) = squeeze(unique_orig_hex(symm,  21:141, :))';
inps(:, :, 2) = squeeze(unique_fitted_hex(symm,  21:141, :))';
inps(find(rel_sess<200), :, :) = nan;
rng(07062023)
ou = run_cluster_based_perm_test(inps, 1, 1000, 1, 97.5, 2);
ou.cluster_stats
plot(ou.cluster_stats.survived_mass_UB - 0.94, 'color', blueCol, 'linewidth', 2)
plot(ou.cluster_stats.survived_mass_LB - 0.94, 'color', redCol, 'linewidth', 2)
legend([a, b], {'Veridical', 'Distorted'}, 'box', 'off')


inps = [];
inps(:, :, 1) = squeeze(unique_orig_hex(3, 1:140, :))';
inps(:, :, 2) = zeros(size(squeeze(unique_orig_hex(3, 1:140, :))'));
inps(find(isnan(unique_orig_hex(3, 1, :))), :, 2) = nan;
ou = run_cluster_based_perm_test(inps, 1, 1000, 1, 95, 2);
ou.cluster_stats


inps = [];
inps(:, :, 1) = squeeze(unique_fitted_hex(3, :, :))';
inps(:, :, 2) = zeros(size(squeeze(unique_fitted_hex(3, :, :))'));
inps(find(isnan(unique_fitted_hex(3, 1, :))), :, 2) = nan;
ou = run_cluster_based_perm_test(inps, 1, 1000, 1, 95, 2);
ou.cluster_stats


subplot(2, 5, 8) % we first obsrved it was significantly stronger in VMPFC in the original time window
rel_rg = [5];
[unique_fitted_hex, unique_orig_hex, unique_hex_diff, rel_sess] = compute_session_level_stats(animal, session, fitted, original,brain_region, rel_rg, nTR);
rel_rws = find(rel_sess > 200);

figInfo = {'VMPFC', 'Grid-like Code', 'Betas [a.u.]', [], [-0.015, 0.1], [1 2], ...
           {'Veridical', 'Distorted'}, [], {}, 14, [], []};
[curr] = plotBar([squeeze(nanmean(unique_orig_hex(3, 51:81, rel_rws), 2)), squeeze(nanmean(unique_fitted_hex(3, 51:81, rel_rws), 2))], figInfo, blueCol, redCol, curr, confInt, 0); % peak neg effect (vmpfc correct)
axis square
[a,b,c,d]=ttest(squeeze(nanmean(unique_orig_hex(:, 51:81, rel_rws), 2))', squeeze(nanmean(unique_fitted_hex(:, 51:81, rel_rws), 2))')
[a,b,c,d]=ttest(squeeze(nanmean(unique_fitted_hex(:, 51:81, rel_rws), 2))')
[a,b,c,d]=ttest(squeeze(nanmean(unique_orig_hex(:, 51:81, rel_rws), 2))')


subplot(2, 5, 9) % we plotted only ACC as it was the only one with a significant distorted signal
rel_rg = [1];
[unique_fitted_hex, unique_orig_hex, unique_hex_diff, rel_sess] = compute_session_level_stats(animal, session, fitted, original,brain_region, rel_rg, nTR);
rel_rws = find(rel_sess > 200);

figInfo = {'ACC', 'Grid-like Code', 'Betas [a.u.]', [], [-0.015, 0.1], [1 2], ...
           {'Veridical', 'Distorted'}, [], {}, 14, [], []};
[curr] = plotBar([squeeze(nanmean(unique_orig_hex(3, 51:81, rel_rws), 2)), squeeze(nanmean(unique_fitted_hex(3, 51:81, rel_rws), 2))], figInfo, blueCol, redCol, curr, confInt, 0); % peak neg effect (vmpfc correct)
axis square
[a,b,c,d]=ttest(squeeze(nanmean(unique_orig_hex(:, 51:81, rel_rws), 2))', squeeze(nanmean(unique_fitted_hex(:, 51:81, rel_rws), 2))')
[a,b,c,d]=ttest(squeeze(nanmean(unique_fitted_hex(:, 51:81, rel_rws), 2))')
[a,b,c,d]=ttest(squeeze(nanmean(unique_orig_hex(:, 51:81, rel_rws), 2))')

subplot(2, 5, 10) % we plotted only ACC as it was the only one with a significant distorted signal
rel_rg = [3];
[unique_fitted_hex, unique_orig_hex, unique_hex_diff, rel_sess] = compute_session_level_stats(animal, session, fitted, original,brain_region, rel_rg, nTR);
rel_rws = find(rel_sess > 200);

figInfo = {'OFC', 'Grid-like Code', 'Betas [a.u.]', [], [-0.015, 0.1], [1 2], ...
           {'Veridical', 'Distorted'}, [], {}, 14, [], []};
[curr] = plotBar([squeeze(nanmean(unique_orig_hex(3, 51:81, rel_rws), 2)), squeeze(nanmean(unique_fitted_hex(3, 51:81, rel_rws), 2))], figInfo, blueCol, redCol, curr, confInt, 0); % peak neg effect (vmpfc correct)
axis square
[a,b,c,d]=ttest(squeeze(nanmean(unique_orig_hex(:, 51:81, rel_rws), 2))', squeeze(nanmean(unique_fitted_hex(:, 51:81, rel_rws), 2))')
[a,b,c,d]=ttest(squeeze(nanmean(unique_fitted_hex(:, 51:81, rel_rws), 2))')
[a,b,c,d]=ttest(squeeze(nanmean(unique_orig_hex(:, 51:81, rel_rws), 2))')



saveImage(curr, ['new_grid_warping', '1'], 'C:\Users\sebas\OneDrive - University College London\Sebastijan\KennerleyLab\Figures\grids_working\paper_figures\revision\', [16000, 10000]);



% grid warping supps 

% comparing acc / ofc at original window 
curr = figure(77)
subplot(1, 3, 1) % ofc
rel_rg = [5];
[unique_fitted_hex, unique_orig_hex, unique_hex_diff, rel_sess] = compute_session_level_stats(animal, session, fitted, original,brain_region, rel_rg, nTR);
rel_rws = find(rel_sess > 200);
plot(1:150, repmat(0, [1, 150]), '--k')
a=plotmse(squeeze(unique_orig_hex(3, 1:end, rel_rws))', blueCol, [0, 1]); hold on
b=plotmse(squeeze(unique_fitted_hex(3, 1:end, rel_rws))', redCol, [0, 1])
figElements(curr, 'VMPFC - Veridical / Distorted Grid-Like Code', 'Time [msec]', 'Betas [a.u.]', [20, 140], [-0.075, 0.1], [24, 41, 81, 121], ...
                {'Cue onset', '-400', 'Choice onset', '+400'}, [], {}, 18, [], []); box off

subplot(1, 3, 2) % acc
rel_rg = [1];
[unique_fitted_hex, unique_orig_hex, unique_hex_diff, rel_sess] = compute_session_level_stats(animal, session, fitted, original,brain_region, rel_rg, nTR);
rel_rws = find(rel_sess > 200);
plot(1:150, repmat(0, [1, 150]), '--k')
a=plotmse(squeeze(unique_orig_hex(3, 1:end, rel_rws))', blueCol, [0, 1]); hold on
b=plotmse(squeeze(unique_fitted_hex(3, 1:end, rel_rws))', redCol, [0, 1])
figElements(curr, 'ACC - Veridical / Distorted Grid-Like Code', 'Time [msec]', 'Betas [a.u.]', [20, 140], [-0.075, 0.1], [24, 41, 81, 121], ...
                {'Cue onset', '-400', 'Choice onset', '+400'}, [], {}, 18, [], []); box off

subplot(1, 3, 3) % acc
rel_rg = [3];
[unique_fitted_hex, unique_orig_hex, unique_hex_diff, rel_sess] = compute_session_level_stats(animal, session, fitted, original,brain_region, rel_rg, nTR);
rel_rws = find(rel_sess > 200);
plot(1:150, repmat(0, [1, 150]), '--k')
a=plotmse(squeeze(unique_orig_hex(3, 1:end, rel_rws))', blueCol, [0, 1]); hold on
b=plotmse(squeeze(unique_fitted_hex(3, 1:end, rel_rws))', redCol, [0, 1])
figElements(curr, 'OFC - Veridical / Distorted Grid-Like Code', 'Time [msec]', 'Betas [a.u.]', [20, 140], [-0.075, 0.1], [24, 41, 81, 121], ...
                {'Cue onset', '-400', 'Choice onset', '+400'}, [], {}, 18, [], []); box off

saveImage(curr, ['new_grid_warping_supp', '1'], 'C:\Users\sebas\OneDrive - University College London\Sebastijan\KennerleyLab\Figures\grids_working\paper_figures\revision\', [12000, 3000]);



end



