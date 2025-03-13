clear; clf;

%% load data structs

% cd(go into data folder) 
rel_files = {'overallRipples', 'auxStruct', 'nnmf_example_chan'};
for fl = 1:length(rel_files)
    load(rel_files{fl});
end

%%
curr = figure(7);

foo  = preprocessLFP(nnmf_example_chan(1).ch, 0, auxStruct.rippleBounds, auxStruct.Fs, 20); 
foo2 = preprocessLFP(nnmf_example_chan(1).ch, 0, auxStruct.widebandBounds_supps, auxStruct.Fs, 20); 
chan = 233;
currChan = zeros(size(overallRipples(chan).trialRipplesHz));
currChanRipple = overallRipples(chan).trialRipples; currChanRipple = currChanRipple(~isnan(currChanRipple(:, 1)), :);
currChanRipple = currChanRipple(currChanRipple(:, 1)~=0, :);
currChanRipple = currChanRipple(currChanRipple(:, 4) > 2000, :);
trials = currChanRipple(:, 1)';
currChanRipple = currChanRipple(:, 5:end);
currChanMorlet                    = preprocessLFP(nnmf_example_chan(1).ch, 1, auxStruct.widebandBounds_supps, auxStruct.Fs, 20);
frequencies                       = currChanMorlet.frequencies;
currChanMorlet                    = currChanMorlet.lfpMagnitude;

trIndx = find(trials == 156);
SNR           = nan(200, length(frequencies), length(trials));
rippleMaxFreq = nan(length(trials), 1);
jj = 1;
pureRipples = nan(300, 301);
relTrials   = nan(300, 1);
origTrials  = nan(size(currChan, 1), 1);
togglePlotting = 1;

tr = 156;
tmpSNR = []; projections = [];

trCurr = currChanRipple(trIndx, (~isnan(currChanRipple(trIndx, :))));

trCurrOrig = trCurr; trCurrOrig = trCurrOrig(trCurrOrig > 0);
ripplePower = squeeze(currChanMorlet(tr, trCurrOrig, :));
[a,b]=nnmf(ripplePower, auxStruct.nClusters, 'replicates', auxStruct.REPLICATES); % b is factor loading
denom = sum(ripplePower.^2);
for k = 1:auxStruct.nClusters
    projections(:, :, k) = a(:, k)*b(k, :);
end
numer = squeeze(sum((projections - ripplePower).^2));
for k = 1:auxStruct.nClusters
    tmpSNR(:, k) = denom./numer(:, k)';
end

factor = max(tmpSNR); [mxFactor, mxId] = max(factor);
bestFactor = tmpSNR(:, mxId);

subplot(3, 5, 4)
args = {1:size(ripplePower, 1), frequencies, ripplePower'};
surf(args{:},'edgecolor','none'); axis square
view(0,90); colormap viridis
axis tight;
figElements(curr, 'Original Power Spectrum', '', 'Frequency', [], [], [], {}, [], {}, 14, 1, []);

subplot(3, 5, 3)
plot(1:200, repmat(5, [200, 1]), ':k');
plot(frequencies, tmpSNR, 'linewidth', 3); box off
legend('Cluster 1', 'Cluster 2', 'Cluster 3'); hold on; 
figElements(curr, '', 'Frequencies', 'Signal to Noise ratio [a.u.]', [], [], [], {}, [], {}, 14, [], []);

subplot(3, 5, 5)
multProj = a(:, mxId)*b(mxId, :);

args = {1:size(ripplePower, 1), frequencies, multProj'};
surf(args{:},'edgecolor','none'); axis square
view(0,90); colormap viridis
axis tight;
figElements(curr, 'Reconstructed Power Spectrum', '', 'Frequency', [], [], [], {}, [], {}, 14, 1, []);

subplot(3, 5, 1)
plot(foo.lfp(tr, :)', 'k'); box off; hold on
plot(trCurrOrig, foo.lfp(tr, trCurrOrig)', 'r'); box off
figElements(curr, '', '', '', ...
    [1 2500], [-0.08, 0.08], [], {}, [], {}, 14, 1, []);
set(gca, 'visible', 'off')

subplot(3, 5, 2)
plot(foo2.lfp(tr, :), 'k'); box off; hold on
plot(trCurrOrig, foo2.lfp(tr, trCurrOrig)', 'r'); box off
figElements(curr, '', '', '', ...
    [1 2500], [-0.1, 0.1], [], {}, [], {}, 14, 1, []);
set(gca, 'visible', 'off')

