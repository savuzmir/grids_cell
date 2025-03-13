function out = figure5_rippleSupplements()
clear; clf;

%% load data structs
% cd(go into data folder) 
rel_files = {'cleanedRipples_Full', 'temporal_firing_supp', 'sample_lfp_chans_supps', 'auxLFP', 'sample_lfp_chans', 'channel_info', 'ripple_aligned_fr', 'auxStruct'};
for fl = 1:length(rel_files)
    load(rel_files{fl});
end

%% generate single channel samples from 1 session

curr = figure(7)
for i = 1:9
    subplot(9, 1, i)
    plot(filtered_out(i, :), 'k'); hold on
    set(gca, 'visible', 'off')
end

%% vmpfc across all channels 

vmpfcAverage = [];
vmpfcIndx = 1;
vmpfc_channels = find(auxStruct.brain_region_channels == 5);
for i = vmpfc_channels
    currChan = cleanedRipples(i).trialRipplesHzPure;
    currChan = currChan(:, 2:end);
    vmpfcAverage(vmpfcIndx, :) = nanmean(currChan);
    vmpfcIndx = vmpfcIndx + 1;
end

curr = figure(6)
subplot(3, 4, 5)
[~, idx] = sort(vmpfcAverage(:, 800));
imagesc(vmpfcAverage(idx , 1:1600)*100)
colormap viridis
caxis([0, 10]); axis square
h = colorbar;
ylabel(h, 'Probability [%]')

figElements(curr, '', '', 'VMPFC Channels', ...
            [], [], [400, 800, 1200], {'-400', 'Choice ON', '400'}, [], {}, 24, [], []);

%% plot firing rate information

% now join them together
allPeak             = [];
allPeakBaseline     = [];
allPeakControl = [];

for i = 1:auxStruct.num_vmpfc_cells

    allPeak            = [allPeak; nanmean(data_file(i).container)];
    allPeakControl     = [allPeakControl; nanmean(data_file(i).container_NR)];
    allPeakBaseline    = [allPeakBaseline; nanmean(data_file(i).baselineFR)];

end

offset = 175; % window centred on ripple event; need righ offset
a      = q_smooth(allPeak(:, (auxStruct.maxIndx-(offset - 1)):(auxStruct.maxIndx + offset)), 20, 1) - allPeakBaseline;
b      = q_smooth(allPeakControl(:, (auxStruct.maxIndx-(offset - 1)):(auxStruct.maxIndx + offset)), 20, 1) - allPeakBaseline;

[~,p, ~, t] = ttest(nanmean(a(:, auxStruct.ripple_fr_sup_roi), 2), nanmean(b(:, auxStruct.ripple_fr_sup_roi), 2))
% 3.36, p = 9.67e-04

subplot(3, 4, 6)
figInfo = {'', '', 'Fr [Hz]', [], [], [1 2], {'Ripple', 'No ripple'}, [], {}, 24, [], []};
plotBar([nanmean(a(:, auxStruct.ripple_fr_sup_roi), 2), nanmean(b(:, auxStruct.ripple_fr_sup_roi), 2)]*1000, figInfo, auxStruct.blueVec, auxStruct.yellowVec, curr, [0, 1], 0);  axis square

h = [];
for k = 3:4 % loop through both
    subplot(3, 4, 4+k)

    for j = 1:4
        if k == 4
            h{j} = plotmse(temporal_firing_supp.frTime_corrected(:, :, j) * 1000, auxStruct.areaCol(j, :), [0, 1]); hold on; axis square
            figElements(curr, '', 'Time [msec]', 'Firing Rate_{\Delta}[Hz]', [1 141], [], [41 81 121], {'-400', 'Choice ON' '+400'}, [], {}, 24, [], [])
        else
            h{j} = plotmse(temporal_firing_supp.frTime(:, :, j) * 1000, auxStruct.areaCol(j, :), [0, 1]); hold on; axis square
            figElements(curr, '', 'Time [msec]', 'Firing Rate [Hz]', [1 141], [], [41 81 121], {'-400', 'Choice ON' '+400'}, [], {}, 24, [], [])
        end
    end
end


end
