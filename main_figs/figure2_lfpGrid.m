function figure2_lfpGrid(auxStruct)

%% load data structs
% cd(go into data folder) 
rel_files = {'auxStruct', 'auxLFP', 'out_results_d1', 'out_results_d2'};
for fl = 1:length(rel_files)
    load(rel_files{fl});
end

%% look at session level stats for dataset 1

% generate stats
for rg = 1:length(auxStruct.areaCodes)
    currRg = auxStruct.relRegs{rg};

    % plot averaged stats
    avgBetasPE  = [];
    avgBetasSR  = [];
    avgSignalPE = [];
    dy = [];

    indx = 1;
    sessions = auxStruct.chanSessions(currRg, :);

    for an = 0:1

        currSessions = sessions(sessions(:, 6) == an, 2);
        tt = auxStruct.dayInformation_stimSet(sessions(:, 1));
        tt = tt(sessions(:, 6) == an);

        unqSessions = unique(currSessions);
        rlReg = intersect(currRg, sessions(sessions(:, 6) == an, 4));

        betasPE   = squeeze(nanmean(out_dat1.RegressionResults_glmfit.Test(1, rlReg, :, :), 4));
        betasSR   = squeeze(nanmean(out_dat1.TestReal(rlReg, :, :, :), 3));
        avgSignal = out_dat1.signal(:, rlReg);

        for us = unqSessions'
            avgBetasPE(indx, :)    = nanmean(betasPE(currSessions == us, :), 1);
            avgBetasSR(indx, :, :) = squeeze(nanmean(betasSR(currSessions == us, :, :), 1));
            avgSignalPE(indx, :)   = squeeze(nanmean(avgSignal(:, currSessions == us), 2)');
            dy(indx)               = median(tt(currSessions == us));
            indx = indx + 1;
        end

    end
    
    % stats for regions
    [a,b,c,d]=ttest(squeeze(avgBetasPE));
    [d.tstat; b; b.*5] % raw + bonf corrected pvals

    sesBetasFullPE{rg} = avgBetasPE;
    sesBetas{rg}       = avgBetasPE(:, :);
    sesSignalPE{rg}    = avgSignalPE(:, :);
    sesBetasSR{rg}     = squeeze(avgBetasSR(:, 3, :));
end

% look at cluster extend significance
[~, ~, ~, trueStats] = ttest(squeeze(avgBetasSR(:, 3, :)));

permStatsSR = [];
for np = 1:(auxStruct.nPerm + 1)
    for rg = 1:4
        [~,~,~,tstat] = ttest(squeeze(out_dat1.avgBetasPermSRFull(:, rg, :, np)));
        permStatsSR(:, rg, np) = tstat.tstat;
    end
end

permStatsRegSR          = nan(1, size(permStatsSR, 1), (auxStruct.nPerm + 1));
permStatsRegSR(1, :, :) = squeeze(permStatsSR(:, 4, :));

[out_vals]=compute_cluster_stats(trueStats.tstat(1, auxStruct.rel_sr_tp), permStatsRegSR(:, auxStruct.rel_sr_tp, :), 95, 3);

%% generate full plot

font_sz = 24;

% main sixfold plot
confInt = [0, 0];
curr=figure(15)
subplot(3, 5, 1)
figInfo = {'VMPFC: sessions', 'Symmetries', 'Betas [a.u.]', [], [-0.06, 0.12], [1 2 3 4 5], {'4', '5', '6', '7', '8'}, [], {}, font_sz, [], []};
[curr] = plotBar(avgBetasPE, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, confInt, 0);

% other regions in same time period
subplot(3, 5, 2:3)
h = [];
onsetModif = 0;
onsetBin = 81 - onsetModif;
endModif   = 141;
nTimepoints = length(1:(endModif-onsetModif));

%plotmse(squeeze(sesBetasSR{rg}), auxStruct.areaCol(rg, :), [0, 1]); hold on
rectangle('Position', [51, -0.1, 30, 0.215], ...
    'FaceColor', [1/255, 114/255, 178/255, 0.2], ...
    'EdgeColor', 'none'); hold on
plot(1:(endModif-onsetModif), repmat(0, [1, nTimepoints]), ':k', 'linewidth', 0.6); hold on
plot(repmat(onsetBin, 2, 1), [-0.1, -0.09], 'linewidth', 5, 'color', 'k');
pvalsOver = [];
for rg = 1:4
    h(rg)=plotmse(squeeze(sesBetasSR{rg}), auxStruct.areaCol(rg, :), [0, 1]); hold on

    [~,b,~,~]=ttest(squeeze(sesBetasSR{rg}));

    fo = nan(1, 161);
    fo(b<0.05) = 0.115 + rg/200;
    plot(1:size(permStatsSR, 1), fo, 'color', auxStruct.areaCol(rg, :), 'linewidth', 2);
end

figElements(curr, '', 'Time [msec]', 'Betas [a.u.]', [1 141], [-0.1, 0.16], [41-onsetModif 81-onsetModif 121-onsetModif], ...
    {'-400', 'Choice ON', '400'}, [], {}, font_sz, [], []); box off

subplot(3, 5, 4:5)
confInt = [0, 0.975];
figInfo = {'VMPFC: sessions', 'Angles', 'z-score [a.u.]', [], [-0.12, 0.15], [1 3 5 7 9 11], {'\omega', '\omega + 60', '\omega + 120', '\omega + 180', '\omega + 240', '\omega + 300'}, [], {}, font_sz, [], []};
[curr] = plotBar(sesSignalPE{4}, figInfo, auxStruct.pinkVec, auxStruct.grayVec, curr, confInt, 0);

% select rel session
i = unique(currSessions)';
i = i(4);

subplot(3, 4, 6:7)
figInfo = {'Session: 46', 'Angles', 'z-score [a.u.]', [], [-0.45, 0.45], [1 3 5 7 9 11], {'\omega', '\omega + 60', '\omega + 120', '\omega + 180', '\omega + 240', '\omega + 300'}, [], {}, font_sz, [], []};

xob = squeeze(avgSignal(:, currSessions == i));
[curr] = plotBar(xob', figInfo, auxStruct.pinkVec, auxStruct.grayVec, curr, confInt, 0);

tmpOn  = squeeze(nanmean(xob(1:2:end, :)));
tmpOff = squeeze(nanmean(xob(2:2:end, :)));
currOnOffJoint = [tmpOn', tmpOff'];

subplot(3, 4, 8)
figInfo = {'Session: 46', 'Angles', 'z-score [a.u.]', [], [-0.25, 0.25], [1 2], {'Aligned', 'Unaligned'}, [], {}, font_sz, [], []};
[curr] = plotBar(currOnOffJoint, figInfo, auxStruct.pinkVec, auxStruct.grayVec, curr, confInt, 0);

subplot(3, 4, 5)
% look at nr significant sessions
uppVMPFC = squeeze(out_dat1.uppThreshold(:, :, 4));
uppVMPFC = uppVMPFC(1:16, :);
percSig = sum(avgBetasPE > uppVMPFC)./size(uppVMPFC, 1);
box off
b = bar(1:5, [percSig]*100, 'edgeColor', 'k', 'linewidth', 3); hold on
for i = 1:length(b)
    b(i).FaceColor = 'flat';

    b(i).CData = auxStruct.greenVec;
end

% sig sess w bonf corr
BinomTest(sum(avgBetasPE > uppVMPFC), size(avgBetasPE, 1), auxStruct.binom_thresh).*5

b.FaceAlpha = 0.7;
figElements(curr, 'Significant Sessions', 'Symmetries', 'Percentage [%]', [], [1, 50], [1 2 3 4 5], ...
    {'4', '5', '6', '7', '8'}, [], {}, font_sz, [], []); box off

%% addon independent dataset

curr = figure(7)
% main sixfold plot
confInt = [0, 0];
subplot(3, 5, 11)
figInfo = {{'Independent Dataset'}, 'Symmetries', 'Betas [a.u.]', [], [-0.06, 0.12], [1 2 3 4 5], {'4', '5', '6', '7', '8'}, [], {}, font_sz, [], []};
[curr] = plotBar(out_dat2.vmpfc_pe, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, confInt, 0);

% sixfold effect in dat 2
[a,b,c,d]=ttest(out_dat2.vmpfc_pe(:, 3))

% aux plot for other regions in same time period
subplot(3, 5, 12:13)
h           = [];
onsetBin    = 151;
endModif    = 291;
nTimepoints = length(1:(endModif));

rectangle('Position', [161, -0.1, 40, 0.215], ...
    'FaceColor', [1/255, 114/255, 178/255, 0.2], ...
    'EdgeColor', 'none'); hold on
plot(1:(endModif), repmat(0, [1, nTimepoints]), ':k', 'linewidth', 0.6); hold on
plot(repmat(onsetBin, 2, 1), [-0.1, -0.09], 'linewidth', 5, 'color', 'k');

b        = plotmse(out_dat2.ofc_sr(:, 1:290, 3), auxStruct.areaCol(3, :), [0, 1]);
a        = plotmse(out_dat2.vmpfc_sr(:, 1:290, 3), auxStruct.areaCol(4, :), [0, 1]);

% add sig traces for visualisation
[~,pvl] = ttest(out_dat2.vmpfc_sr(:, :, 3));
rel_a = nan(290, 1);
rel_a(find(pvl<.05)) = 0.1;
[~,pvl] = ttest(out_dat2.ofc_sr(:, :, 3));
rel_b = nan(290, 1);
rel_b(find(pvl<.05)) = 0.08;

plot(1:290, rel_a, 'linewidth', 2,  'color', auxStruct.areaCol(4, :));
plot(1:290, rel_b, 'linewidth', 2,  'color', auxStruct.areaCol(3, :));

figElements(curr, 'Independent Dataset', 'Time [msec]', 'Betas [a.u.]', [30, 291], [], [70, 150, 220], {'Option 1 ON', 'Option 2 ON', 'Choice ON'}, [], {}, font_sz, [], []);
legend([b, a], {'OFC', 'VMPFC'}); legend box off

curr_on_off = out_dat2.curr_on_off;
confInt = [0, 0.975];

subplot(3, 5, 14:15)
figInfo = {{'', 'Independent Dataset'}, 'Angles', 'z-score [a.u.]', [], [-0.15, 0.15], [1 3 5 7 9 11], {'\omega', '\omega + 60', '\omega + 120', '\omega + 180', '\omega + 240', '\omega + 300'}, [], {}, font_sz, [], []};
[curr] = plotBar(curr_on_off, figInfo, auxStruct.pinkVec, auxStruct.grayVec, curr, confInt, 0);

% just on-off; very significant in t-test
ons  = curr_on_off(:, 1:2:12);
offs = curr_on_off(:, 2:2:12);
rel = [reshape(ons, [6*13, 1]), reshape(offs, [6*13, 1])];
[a,b,c,d]=ttest(rel(:, 1), rel(:, 2))


end