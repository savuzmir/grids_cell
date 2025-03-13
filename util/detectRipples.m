function overallRipples = detectRipples(inp, currLFP, tot_len, lfpAverage, filter_tp)

% detects ripples

% inp is settings, currLFP is LFP data, auxLFP is aux data

overallRipples = struct;
overallRipples.trialRipples = [];

% get gaussian smoothing kernel
gaussFilter = auxStruct.gaussFilter;

cutoff       = inp.sdev;
modif        = inp.minDiff;
windLengthLB = inp.lb;
windLengthUB = inp.ub;
cutoffLength = inp.cutoffLength;

rippleBounds = inp.rippleBounds;

tic

% bandpass in ripple band
rippleFilt = preprocessLFP_ripple(currLFP(1).ch, 0, [rippleBounds(1), rippleBounds(2)], inp.Fs, filter_tp);

% define a container for the ripples
overallRipples.trialRipples = nan(350, 350); % this assumes a max of 350 evs x 350 msec which is a reasonable bound for here 

trueRippleIndx = 1;

% define nr of trials
nTrials = size(rippleFilt.lfp, 1);

% get global session estimates of mean and stdev by concatenating all trials
sessionConcat = [];
% we define the global mean and standard deviation based on concatenated epoched data
for tr = 1:nTrials
    sessionConcat = [sessionConcat, rippleFilt.lfp(tr, :)];
end

% compute envelope
envelope = abs(hilbert(sessionConcat));

% smooth it with the filter
smoothedEnvelope = conv(envelope,gaussFilter,'same');

if filter_tp == 1
    % take global mean and std to keep consistent across comparisons
    men = lfpAverage.mn; 
    stdev = lfpAverage.std;
else
    men   = nanmean(smoothedEnvelope);
    stdev = nanstd(smoothedEnvelope);
end

% now we use the individual trials estimates
absIndx = 1;
for tr = 1:nTrials

    % determine envelope
    envelope = abs(hilbert(rippleFilt.lfp(tr, :)));

    % smooth
    smoothedEnvelope = conv(envelope,gaussFilter, 'same');

    % zscore the envelope wrt global values
    smoothEnvZ = (smoothedEnvelope - men)/stdev;

    % find all events based on a preliminary threshold pass like in Karlsson & Frank
    % from these needs to exceed cuttoff and length criteria defined below
    candidateEvents = regionprops(smoothEnvZ >= 0, 'Area', 'PixelIdxList');

    % find candidate ripples; we predefine an arbitrary container
    trueCandidates = nan(20, 500);
    concatEvents = [];

    for ce = 1:size(candidateEvents, 1)
        % check their properties
        timeIndices = [candidateEvents(ce).PixelIdxList];
        currCand = smoothEnvZ(timeIndices);

        % find peak of currCand
        [~, mx] = max(currCand);

        % put back into correct time frame

        refAligned = timeIndices(mx);

        % take lb length / 2 on each side and create the adjusted bounds
        lbAdjusted = floor(refAligned - (cutoffLength - 1)/2);
        ubAdjusted = ceil(refAligned + (cutoffLength - 1)/2);

        peakBounds = lbAdjusted:ubAdjusted;

        % when this evaluates to true, they are at the edges and we cannot evaluate them properly
        if any(peakBounds <= 0)  || any(peakBounds > tot_len)
            ...
        else

            % current cand has to be > cutoff, and its length needs to pass the criteria
            currCandPeak = smoothEnvZ(peakBounds);
            if  all(currCandPeak>=cutoff) & ((length(currCand) >= windLengthLB) & (length(currCand) <= windLengthUB)) % upper bound removes clipping at edges
                concatEvents = [concatEvents candidateEvents(ce).PixelIdxList'];
            end

        end
    end

    % merge close peaks
    peakDiffs = diff(concatEvents);
    peakDiffs(diff(concatEvents) < modif) = 1;
    peakDiffs = find(peakDiffs > 1); % check how many remain

    if ~isempty(peakDiffs)
        for trCe = 1:(size(peakDiffs, 2) + 1)

            if trCe == 1 % if its the first
                rippleLength = concatEvents(1:peakDiffs(trCe));
            elseif trCe == 2 & size(peakDiffs, 2) == 2 % if there are 3 ripples
                rippleLength = concatEvents((peakDiffs(trCe-1) + 1):peakDiffs(trCe));
            elseif trCe == 2 & size(peakDiffs, 2) == 1 % if there are 2 ripples then go to the end
                rippleLength = concatEvents((peakDiffs(trCe-1) + 1):length(concatEvents));
            elseif trCe == (size(peakDiffs, 2) + 1) % for final ones; note this may not work on very longer windows
                rippleLength = concatEvents((peakDiffs(trCe-1) + 1):length(concatEvents));
            else
                rippleLength = concatEvents((peakDiffs(trCe-1) + 1):peakDiffs(trCe));
            end

            trueCandidates(trCe, 1:length(rippleLength)) = rippleLength;
        end
    else
        trueCandidates(1, 1:length(concatEvents)) = concatEvents;
    end

    % remove nans
    trueRows = sum(~isnan(trueCandidates(:, 1)));

    % how many discrete events do we have
    rippleContainer = zeros(50, tot_len);
    for rw = 1:trueRows

        currRipple = trueCandidates(rw, :);
        currRipple = currRipple(~isnan(currRipple));

        currRipple = currRipple(currRipple ~=0);

        % store time index of max peak
        [~,timeInd]= max(smoothEnvZ(currRipple));
        overallRipples.trialRipples(absIndx, 1:(length(currRipple) + 2)) = [tr, currRipple(timeInd), currRipple];

        rippleContainer(rw, currRipple) = 1;
        absIndx = absIndx + 1;
    end

    overallRipples.trialRipplesHz(tr, 1:(tot_len + 1)) = [tr, sum(rippleContainer)];

end

end