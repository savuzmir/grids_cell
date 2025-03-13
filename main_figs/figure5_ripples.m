function [curr, newStats] = figure5_ripples(ripplesFixation, overallRipples, auxLFP, auxStruct)
clear; clf;

% cd(go into data folder) 
rel_files = {'sample_lfp_chans', 'auxStruct', 'AuxiliaryCells', 'allUnits', 'auxLFP', 'sample_lfp_chans', 'channel_info', 'ripple_aligned_fr'};
for fl = 1:length(rel_files)
    load(rel_files{fl});
end

%% plot pipeline examples for ripples

% start figure
curr=figure(3);

% rel examples
full_fp_dlpfc_a = preprocessLFP(rel_lfp_chans(1).ch, 0, [1, 250], auxStruct.Fs, 14);
full_fp_dlpfc_b = preprocessLFP(rel_lfp_chans(2).ch, 0, [1, 250], auxStruct.Fs, 14);
full_fp_vmpfc_a = preprocessLFP(rel_lfp_chans(3).ch, 0, [1, 250], auxStruct.Fs, 14);
full_fp_vmpfc_b = preprocessLFP(rel_lfp_chans(4).ch, 0, [1, 250], auxStruct.Fs, 14);
ripp_fp_vmpfc_b = preprocessLFP(rel_lfp_chans(4).ch, 0, [80, 180], auxStruct.Fs, 14); 

j = 177; % trial example

subplot(8, 3, 1:2)
plot(full_fp_vmpfc_a.lfp(j, :), 'k', 'linewidth', 1.5)
figElements(curr, '', '', '', ...
[], [-0.175, 0.175], [], {}, [], {}, 24, [], []);
% VMPFC - ch146

axis off 
subplot(8, 3, 4:5)
plot(full_fp_vmpfc_b.lfp(j, :), 'k', 'linewidth', 1.5); hold on
figElements(curr, '', '', '', ...
[], [-0.175, 0.175], [], {}, [], {}, 24, [], []);
axis off 
% VMPFC - ch147

set(gca, 'visible', 'off')

subplot(8, 3, 13:14)
plot(full_fp_dlpfc_a.lfp(j, :), 'k', 'linewidth', 1.5)
figElements(curr, '', '', '', ...
[], [-0.175, 0.175], [], {}, [], {}, 24, [], []);
axis off 
% DLPFC - ch137

subplot(8, 3, 16:17)
plot(full_fp_dlpfc_b.lfp(j, :), 'k', 'linewidth', 1.5)
figElements(curr, '', '', '', ...
[], [-0.175, 0.175], [], {}, [], {}, 24, [], []);
axis off 

% DLPFC - ch143

% determine envelope
envelope = abs(hilbert(ripp_fp_vmpfc_b.lfp(j, :)));
% we smooth the envelope now

smoothedEnvelope = conv(envelope, auxStruct.gaussFilter,'same');
% zscore the envelope
smoothEnvZ = zscore(smoothedEnvelope);

subplot(8, 3, 10:11)
plot(1:2500, smoothEnvZ(1:2500), 'k', 'linewidth', 1.5); box off, hold on
plot(1:2500, repmat(3, [1, 2500]), '--', 'color', [155, 155, 155]./255)
figElements(curr, '', '', '', ...
[], [], [], {}, [], {}, 24, [], []);
axis off 

% Smoothed Envelope - ch147
subplot(8, 3, 7:8)
plot(ripp_fp_vmpfc_b.lfp(j, :)', 'k', 'linewidth', 1.5); box off
figElements(curr, '', '', '', ...
[], [], [], {}, [], {}, 24, [], []);
axis off 

curr = figure(67)
% Ripple Band - ch147

lowpassAvg  = channel_info.channelAvg(:, (auxStruct.maxIndx-(auxStruct.offset - 1)):(auxStruct.maxIndx + auxStruct.offset));
widebandAvg = channel_info.channelAvg_full(:, (auxStruct.maxIndx-(auxStruct.offset - 1)):(auxStruct.maxIndx + auxStruct.offset));   

subplot(5, 3, 3)
plotmse(lowpassAvg, 'k', [0, 1]); hold on; axis off
axis square 
figElements(curr, '', '', 'mV', ... 
            [126 475], [-0.03, 0.03], [200, 250, 300, 350, 400], {'-100', '-50', 'Peak', '50', '100'}, [], {}, 24, [], []);   

curr = figure(6)
plotmse(widebandAvg, 'k', [0, 1]); hold on
figElements(curr, '', 'Time [msec]', 'mV', [126 475], [-0.03, 0.03], [200, 250, 300, 350, 400], {'-100', '-50', 'Peak', '50', '100'}, [], {}, 24, [], []);   

% now join them together
allPeak             = [];
allPeakBaseline     = [];

num_vmpfc_cells = length(find(auxStruct.brain_region_cells == 5));
for i = 1:num_vmpfc_cells
    allPeak            = [allPeak; nanmean(data_file(i).container, 1)];
    allPeakBaseline    = [allPeakBaseline; nanmean(data_file(i).baselineFR)];
end

curr = figure(67)
offset = 175; % need different cutoff due to different saving 
smoothPeak = q_smooth(allPeak, 20, 1); smoothPeak = smoothPeak(:, (auxStruct.maxIndx-(offset - 1)):(auxStruct.maxIndx + offset)); % need slighlty different to clip better
plot([1, 350], [0, 0], ':k', 'linewidth', 0.5); hold on
plotmse(smoothPeak*1000 - allPeakBaseline*1000, auxStruct.blueVec, [0, 1]); hold on; box off
figElements(curr, '', '', 'Firing Rate [Hz]', [1, 350], [-2, 4.5], [1, 50, 100, 150, 200], {'-100', '-50', 'Peak', '50', '100'}, [], {}, 24, 1, []);   
sz = 4000;


foo = preprocessLFP(rel_lfp_chans(1).ch, 1, [2, 200], 1000, 20); % octaves was set to 20 on old pc 
frequencies = foo.frequencies;

curr = figure(7)
args = {1:300, frequencies', squeeze(nanmean(channel_info.channelAvg_spectral(:, 1051:1350, :)))'};
surf(args{:},'edgecolor','none');
view(0,90);
h = colorbar;
h.Location = 'northoutside';
ylabel(h, 'Power [a.u.]')
caxis([0.0002, 0.00035])
colormap viridis
figElements(curr, '', '', 'Frequency', [1, 300], [], [50, 100 150, 200, 250], {'-100', '-50', 'Ripple Peak', '50', '100'}, [], {}, 24, 1, []);


%% compute relationship between ripple frequency and behaviour

resp = load('cleanedRipples_Full.mat');
rew  = load('cleanedRipplesReward_Full.mat');
fix  = load('cleanedRipplesFixation_Full.mat');

resp = resp.cleanedRipples;
rew  = rew.cleanedRipples;
fix  = fix.cleanedRipples;

for i = 1:length(auxStruct.relRegs)

    currCol = auxStruct.areaCol(i, :);
    rlReg = find(auxStruct.brain_region_channels == auxStruct.areaCodes(i)); 

    % first average
    indJ = 1;
    respAvg = []; fixAvg = []; rewAvg = [];

    for j = rlReg
        respAvg(indJ, :) = nanmean(resp(j).trialRipplesHzPure(:, :))*100;
        fixAvg(indJ, :)  = nanmean(fix(j).trialRipplesHzPure(:, :))*100;
        rewAvg(indJ, :)  = nanmean(rew(j).trialRipplesHzPure(:, :))*100;
        indJ = indJ + 1;
    end
    chanAverages(i).respAvg = respAvg;
    chanAverages(i).fixAvg = fixAvg;
    chanAverages(i).rewAvg = rewAvg;

end

curr = figure(7)
% ripple frequency
h = [];

for i = 1:4
    currCol = auxStruct.areaCol(i, :);

    subplot(2, 3, 4)
    plotmse(chanAverages(i).fixAvg(:, 50:500), currCol, [0, 1]); hold on
    figElements(curr, '', '', 'Mean Ch. Norm. Probability [%]', ...
        [1, 451], [0, 2.5], [50, 450], {'-400', 'Cue ON'}, [], {}, 26, [], []);
    title('Fixation')
    subplot(2, 3, 5)
    plotmse(chanAverages(i).respAvg(:, 200:1200), currCol, [0, 1]); hold on
    figElements(curr, '', '', '', ...
        [], [0, 5], [200, 600], {'-400', 'Choice ON'}, [], {}, 26, [], []);
    title('Choice')
    subplot(2, 3, 6)
    h(i) = plotmse(chanAverages(i).rewAvg(:, 800:1950), currCol, [0, 1]); hold on
    figElements(curr, '', '', '', ...
        [], [0, 2.5], [200 1000], {'Reward ON', '800'}, [], {}, 26, [], []);
    title('Reward')
end

legend(h(:), {'ACC', 'DLPFC', 'OFC', 'VMPFC'}, 'NumColumns', 2, 'location', 'best'); legend box off
sz = 16500;

%% do perm tests for behaviour signals

rng(29042020)

newStats = struct;
relReg = find(auxStruct.brain_region_channels == 5);

cleanedRipples = fix;

for i = relReg
    ripplesFixation(i).trialRipplesHz = cleanedRipples(i).trialRipplesHzPure;
end


% does rewarded on prev trial influence ripple prob in fixation

% split by half
relIndx = 1;
avgRipples = [];
for i = relReg
    avgRipples(relIndx, :) = nanmean(ripplesFixation(i).trialRipplesHz);
    relIndx = relIndx + 1;
end

avgRipples = avgRipples(:, 4:1500);

% take average within fixation and chans with above avg freq
fixAvg = nanmean(avgRipples(:, 1:500), 2);
gr1    = find(fixAvg > median(fixAvg));
relReg = relReg(gr1); 

% need to take clean windows
timing1 = 4:1003;
timing2 = 21:521;

averageRipple = nan(length(relReg), length(timing1), 2, auxStruct.nPerm + 1);
nuissanceInfo = nan(length(relReg), 2);

for np = 1:(auxStruct.nPerm + 1)
    chIndx = 1;
    tic
    for i = relReg

        tmp = ripplesFixation(i).trialRipplesHz(:, timing1);

        if np == 1
            trialDivider = circshift(auxLFP(i).ch.Rewarded==1, 1); 
            nuissanceReg = circshift(ripplesFixation(i).trialRipplesHz(:, 3), 1); %need to take 3rd elem (aux info)

            nuissanceInfo(chIndx, 1)        = nanmean(nuissanceReg(trialDivider));
            nuissanceInfo(chIndx, 2)        = nanmean(nuissanceReg(~trialDivider));

        else
            trialDivider = circshift(auxLFP(i).ch.Rewarded==1, 1);
            trialDivider = shuffle(trialDivider);
        end

        averageRipple(chIndx, :, 1, np) = nanmean(tmp(trialDivider, :));
        averageRipple(chIndx, :, 2, np) = nanmean(tmp(~trialDivider, :));

        chIndx = chIndx + 1;
    end
    toc
end

permStats = [];

for i = 2:(auxStruct.nPerm + 1)
    [~,~,~,tt] = ttest(averageRipple(:, timing2, 1, i), averageRipple(:, timing2, 2, i));
    permStats(:, :, i) = tt.tstat;
end

[~,~,~,stats]=ttest(averageRipple(:, timing2, 1, 1), averageRipple(:, timing2, 2, 1));

out = compute_cluster_stats(stats.tstat, permStats(:, :, 2:end), 97.5, 3);

newStats.rewarded.tstat = stats.tstat;
newStats.rewarded.permTest = out;

%% compute relationship betewen ripple proportion and accuracy

for i = 1:auxStruct.num_channels

    %compute RT for each trial
    currChan = auxLFP(i).ch.Codes;
    RTs = [];
    delib_period = [];

    for tr = 1:size(currChan, 1)

        codeInfo = currChan{tr};
        RTs(tr) = codeInfo(codeInfo(:, 1) == 4, 2) - codeInfo(codeInfo(:, 1) == 62, 2);
        delib_period(tr) = codeInfo(codeInfo(:, 1) == 72, 2) - codeInfo(codeInfo(:, 1) == 4, 2);
    end

    all_RT{i}         = RTs;
    overall_RT(i)     = median(RTs);
    overall_delib(i)  = median(delib_period);

    % define brainer trial
    pay_diff  = sign(auxLFP(i).ch.Left_pay  - auxLFP(i).ch.Right_pay);
    prob_diff = sign(auxLFP(i).ch.Left_prob - auxLFP(i).ch.Right_prob);

    brainer   = (sum(abs([prob_diff pay_diff]), 2) > 1) & (pay_diff ~= prob_diff);

    overall_RT_SPLIT(1, i) = median(RTs(brainer==1));
    overall_RT_SPLIT(2, i) = median(RTs(brainer==0));

    % compute accuracy
    choseLeft = auxLFP(i).ch.ChML > 0;
    EV = (auxLFP(i).ch.Left_pay.*auxLFP(i).ch.Left_prob) > (auxLFP(i).ch.Right_pay.*auxLFP(i).ch.Right_prob);

    accuracy(i) = sum(EV==choseLeft)/length(EV);
end

%% define the quantities for GLM

relReg = find(auxStruct.brain_region_channels == 5); 

% get accuracy
relAccuracy = accuracy(relReg);      relAccuracy = zscore(relAccuracy)';
relRT       = zscore(overall_RT(relReg));

% get subject
relSubject = auxStruct.chanSessions(relReg, 5);

% get day
relDays = auxStruct.dayInformation_stimSet(auxStruct.chanSessions(:, 1), 1); relDays = relDays(relReg);
relSess = auxStruct.chanSessions(:, 2); 
relSess = relSess(relReg);

for i = relReg
    overallRipples(i).trialRipplesHz = resp(i).trialRipplesHzPure;
end

concatTimeRipple = [];
indxR = 1;
for i = relReg
    concatTimeRipple(indxR, :) = nanmean(overallRipples(i).trialRipplesHz(:, 4:end));
    indxR = indxR + 1;
end

wholeRipple = zscore(concatTimeRipple);

overallStats = nan(1, 1200, auxStruct.nPerm);

parfor np = 1:(auxStruct.nPerm + 1)

    if np == 1
        tmpRipple = wholeRipple;
    else
        tmpRipple = wholeRipple(shuffle(1:size(wholeRipple, 1)), :);
    end

    tic
    for i = 1:1200

        currRipple = tmpRipple(:, i+199);

        [~,~,stats] = glmfit([currRipple, relSubject, relRT', relDays], relAccuracy);

        outStats = [stats.t(2)];

        overallStats(:, i, np) = outStats;
    end
    toc
end

trueStats = overallStats(1, :, 1);
permStats = overallStats(1, :, 2:(auxStruct.nPerm + 1));
outOverall = compute_cluster_stats(trueStats, permStats, 97.5, 3)

%      clustLenUB: 129.6250
%      clustLenLB: 164.1500
%     clustMassUB: 321.1497
%     clustMassLB: 455.0636
%          trueUB: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 … ]
%          trueLB: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 … ]
%       trueLenUB: [94 6 203 367]
%       trueLenLB: []
%      trueMassUB: [222.8110 11.7295 792.9433 1.4359e+03]
%      trueMassLB: []
%      trueIndxUB: {[94×1 double]  [6×1 double]  [203×1 double]  [367×1 double]}
%      trueIndxLB: {}

%        threshUB: [1.9306 1.9122 1.8991 1.8602 1.9233 1.9348 1.9190 1.9200 1.9473 1.8772 1.8595 1.8857 1.8856 1.8971 1.8811 1.8737 1.9237 1.9171 1.9202 … ]
%        threshLB: [-1.9334 -1.8723 -1.8170 -1.7848 -1.7598 -1.7347 -1.7502 -1.8354 -1.8298 -1.8564 -1.8585 -1.8475 -1.8538 -1.8884 -1.9028 -1.8673 … ]

curr = figure(6)
subplot(1, 3, 1)
currRipple = wholeRipple(:, 457);

[~,~, residX] = glmfit([relSubject, relRT', relDays], currRipple);
[~,~, residY] = glmfit([relSubject, relRT', relDays], relAccuracy);

scatterfit(residX.resid, residY.resid, 0); lsline
axis square
figElements(curr, 'r = .47, p < .001', 'Ripple Probability [a.u.]', 'Accuracy [a.u.]', [], [-1.8, 1.4], [], {}, [], {}, 26, [], []);

curr = figure(67)
subplot(1, 3, 2:3)
h1=plotmse(averageRipple(:, timing2, 1, 1)*100, [0.2738    0.0315    0.3589], [0, 1]); hold on % rewarded
h2=plotmse(averageRipple(:, timing2, 2, 1)*100, [0.7099    0.8688    0.1693], [0, 1])
figElements(curr, '', '', 'Mean Ch. Norm. Probability [%]', [1 501], [0 3], [21 271 481], {'-450', '-250', 'Cue ON'}, [], {}, 26, [], []);
plot(out.trueIndxUB{2}, repmat(2.5, [length(out.trueIndxUB{2}), 1]), 'k', 'linewidth', 4)
legend([h1,h2], 'Rewarded_{trial - 1}', 'Non-rewarded_{trial - 1}', 'location', 'northwest'); legend box off

curr = figure(7)
% plot sig line 
plot([0, 1400], [0, 0], ':k', 'linewidth', 0.6); hold on
plot(trueStats, 'color', auxStruct.greenVec, 'linewidth', 2)
figElements(curr, '', 'Time [msec]', 'T-value [a.u.]', [], [], [200 600, 1000], {'-400', 'Choice ON', '400'}, [], {}, 26, [], []);
plot_sigline_perm(trueStats, outOverall, 6, 1, [0, 1])
%plot(outOverall.threshUB, 'k--', 'linewidth', 0.6); hold on; box off
box off


end