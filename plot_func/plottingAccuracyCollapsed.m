
nBootstraps = 1:size(decoderResults, 1);
timepoints  = 1:size(decoderResults, 2);

% arbitrary choice for plotting
auxStruct.epochDur = auxLFP(1).ch.pre_cue + auxLFP(1).ch.post_cue + 1;
auxStruct.preEvent = auxLFP(1).ch.pre_cue;
auxStruct.windowSize = 100;
auxStruct.slideWidth = 10;
auxStruct.slSizeStart = -400;
auxStruct.slSizeEnd = 800;

bins = [1:auxStruct.slideWidth:(auxStruct.epochDur - auxStruct.slideWidth + 1)];

startBin = round((auxStruct.preEvent + auxStruct.slSizeStart)/auxStruct.slideWidth + 1);
endBin = round((auxStruct.preEvent + auxStruct.slSizeEnd)/auxStruct.slideWidth + 1);

onsetBin = round((auxStruct.preEvent + onset)/auxStruct.slideWidth + 1);


acc = []; corr = []; postProb = []; nTrials = [];
cv = 0;
pp = 0;

for an = 1:2
    for nb = nBootstraps
        for nTP = timepoints
            
            tmp = decoderResults(nb, nTP, an).acc;
            %             tmp = nanmean([tmp{:}]);
            acc(an, nb, nTP) = tmp;
            
            tmp = decoderResults(nb, nTP, an).corr;
            %             tmp = nanmean([tmp{:}]);
            corr(an, nb, nTP) = tmp;
            
            if pp
                for ch = 1:5     % labs, trials, timepoints, ch, an
                    postProbs = nan(5, 80, size(timepoints, 2), 5, 2);
                    
                    if cv
                        for cv = 1:5
                            tmpProb = decoderResults(nb, nTP, an).postProb;
                            tmpProb = tmpProb{cv};
                            tmpCh = decoderResults(nb, nTP, an).chOpt;
                            indxTR = find(tmpCh{cv} == ch);
                            nTR = size(tmpProb(indxTR, :), 1);
                            postProbs(indxTR, 1:5, cv) = tmpProb(indxTR, :);
                        end
                        postProbs = nanmean(nanmean(postProbs), 3);
                        try
                            postProb(nTP, :, an, ch, nb) = postProbs(1, 1:5);
                        catch
                            postProb(nTP, :, an, ch, nb) = nan(1, 5);
                            
                            nTrials(nTP, an, ch, nb) = nTR;
                        end
                    else
                        tmpProb = decoderResults(nb, nTP, an).postProb;
                        tmpCh = decoderResults(nb, nTP, an).chOpt;
                        indxTR = find(tmpCh == ch);
                        nTR = size(tmpProb(indxTR, :), 1);
                        postProbs(1:5, indxTR, nTP, ch, an) = tmpProb(indxTR, :)';
                        
                        try
                            postProb(nTP, :, an, ch, nb) = nanmean(postProbs(1:5, indxTR, nTP, ch, an), 2);
                        catch
                            postProb(nTP, :, an, ch, nb) = nan(1, 5)';
                            
                            nTrials(nTP, an, ch, nb) = nTR;
                        end
                        
                    end
                    
                end
            end
        end
    end
end

%% acc
grayVec = [153/255, 153/255, 153/255];
orangeVec = [230/255, 159/255, 1/255];
blueVec  = [86/255, 180/255, 233/255];


name = 'Prob to Mag VMPFC';

% separate
curr0_a = figure('units','normalized','outerposition',[0 0 1 1]);
plotmse(1:191, squeeze(acc(1, :, :)), grayVec, grayVec, [1, 5], 0.5, 2);
plotmse(1:191, squeeze(acc(2, :, :)), orangeVec, orangeVec, [1, 5], 0.5, 2);
plot(1:nTP, repmat(0.2, [nTP, 1]), '-.k')
ylim([0.10, 0.5])
plot(repmat(onsetBin, [1, 2]), [0, 0.125], 'k', 'linewidth', 3)
figElements(curr0_a, sprintf('Accuracy: %s', name), 'Time [msec]', 'M_{SE} accuracy [%]', [1 nTP], [], [startBin, onsetBin, endBin, endBin + (endBin - onsetBin)], {auxStruct.slSizeStart, 0, auxStruct.slSizeEnd, auxStruct.slSizeEnd*2}, [], {[]}, 12, [], []); hold on

cd('C:\Users\sveselic\OneDrive - University College London\Sebastijan\KennerleyLab\Figures\DecodingNavigation\cue12_collapsed\FR\PROB')
saveas(curr0_a, [name, '.png'])


% joint 
curr0_a = figure('units','normalized','outerposition',[0 0 1 1]);

plotmse(1:191, squeeze(nanmean(acc(:, :, :), 1)), blueVec, blueVec, [1, 5], 0.5, 2);
plot(1:nTP, repmat(0.2, [nTP, 1]), '-.k')
ylim([0.10, 0.5])
plot(repmat(onsetBin, [1, 2]), [0, 0.125], 'k', 'linewidth', 3)
figElements(curr0_a, sprintf('Accuracy: %s', name), 'Time [msec]', 'M_{SE} accuracy [%]', [1 nTP], [], [startBin, onsetBin, endBin, endBin + (endBin - onsetBin)], {auxStruct.slSizeStart, 0, auxStruct.slSizeEnd, auxStruct.slSizeEnd*2}, [], {[]}, 12, [], []); hold on

saveas(curr0_a, [name, '_joint.png'])


    
%% corr
curr0_b = figure(3)
plotmse(1:191, squeeze(corr(1, :, :)), grayVec, grayVec, [1, 5], 0.3, 2);
plotmse(1:191, squeeze(corr(2, :, :)), orangeVec, orangeVec, [1, 5], 0.3, 2);
plot(1:nTP, repmat(0, [nTP, 1]), '-.k')
ylim([-0.2, 0.45])
plot(repmat(onsetBin, [1, 2]), [-0.2, -0.125], 'k', 'linewidth', 3)
figElements(curr0_b, 'Accuracy prob->prob: ACC', 'Time [msec]', 'M_{SE} accuracy [%]', [1 nTP], [], [startBin, onsetBin, endBin, endBin + (endBin - onsetBin)], {auxStruct.slSizeStart, 0, auxStruct.slSizeEnd, auxStruct.slSizeEnd*2}, [], {[]}, 12, [], []); hold on

%% post prob 
postProb = nanmean(postProb, 5); % nTP, :, an, ch, nb
postProbF = squeeze(postProb(:, :, 1, :, :));
postProbM = squeeze(postProb(:, :, 2, :, :));


postProbF = permute(postProbF, [4, 1, 2, 3]);
postProbM = permute(postProbM, [4, 1, 2, 3]);

CUE_RANK = 1;

figure(3)
plotmse(1:191, squeeze(postProbF(:, :, :, CUE_RANK)), blueVec, blueVec, [1, 5], 0.5, 2); legend('Rank 1', 'Rank 2', 'Rank 3', 'Rank 4', 'Rank 5') 
figure(4)
plotmse(1:191, squeeze(postProbM(:, :, :, CUE_RANK)), blueVec, blueVec, [1, 5], 0.5, 2); legend('Rank 1', 'Rank 2', 'Rank 3', 'Rank 4', 'Rank 5') 
