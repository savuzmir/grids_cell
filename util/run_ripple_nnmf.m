function out = run_ripple_nnmf(currLFP, overallRipples, auxStruct, nameType)

for chan = 1:auxStruct.num_channels
    
    tic
    
    currChan = zeros(size(overallRipples(chan).trialRipplesHz));
    currChanGamma = currChan; currChanSigma = currChan;

    % first remove all nan rows to reduce size by checking where first col has nans - this will be true for all of them
    currChanRipple = overallRipples(chan).trialRipples; currChanRipple = currChanRipple(~isnan(currChanRipple(:, 1)), :);
    
    % remove zero rows
    currChanRipple = currChanRipple(currChanRipple(:, 1)~=0, :);
  
    % check how many trials there were
    trials = currChanRipple(:, 1)';
    
    % select only the actual data; note the old convention had additional info in 1st cpl of rows
    ccrf           = currChanRipple;
    currChanRipple = currChanRipple(:, 5:end);

    % get frequency space info
    currChanMorlet                    = preprocessLFP_ripple(currLFP(chan).ch, 1, [2, 200], 1000, 1); 
    frequencies                       = currChanMorlet.frequencies;
    currChanMorlet                    = currChanMorlet.lfpMagnitude;
    
    trIndx = 1;
    SNR           = nan(200, length(frequencies), length(trials));
    rippleMaxFreq = nan(length(trials), 1);
    jj = 1;
    pureRipples = nan(300, 301);
    relTrials   = nan(300, 1);
    origTrials  = nan(size(currChan, 1), 1);

    % for trials where there was a recorded ripple
    for tr = trials 

        tmpSNR = []; projections = [];
        
        % take current trial and from the current trial remove all nans
        trCurr = currChanRipple(trIndx, (~isnan(currChanRipple(trIndx, :))));
        
        % remove also all zeroes because this means they were added in at some point - this shouldnt really happen but is here as a safety check
        trCurrOrig = trCurr; trCurrOrig = trCurrOrig(trCurrOrig > 0);
        
        % take ripple power - we take tr because some ripples occured on same trials with different indices
        ripplePower = squeeze(currChanMorlet(tr, trCurrOrig, :));
        
        [a,b]=nnmf(ripplePower, auxStruct.nClusters, 'replicates', auxStruct.REPLICATES); % b is factor loading

        denom = sum(ripplePower.^2);

        for k = 1:nClusters
            projections(:, :, k) = a(:, k)*b(k, :);
        end

        numer = squeeze(sum((projections - ripplePower).^2));

        % which has highest SNR
        for k = 1:auxStruct.nClusters
            tmpSNR(:, k) = denom./numer(:, k)';
        end
        
        % integrate over factors
        factor = max(tmpSNR); [mxFactor, mxId] = max(factor);

        bestFactor = tmpSNR(:, mxId);
        
        if mxFactor > auxStruct.SNR_THRESH

            [val, idx] = findpeaks(bestFactor, 'MinPeakDistance', 5);

              idx(val < auxStruct.SNR_THRESH) = nan;
             % now remove everything below relevant SNR
              val(val < auxStruct.SNR_THRESH) = nan;
              freqPeaks = frequencies(idx(~isnan(idx)));                 

            if ~(all(freqPeaks >= auxStruct.logo_bounds(1) & freqPeaks < auxStruct.logo_bounds(2)) || all(freqPeaks > auxStruct.logo_bounds(3)) || all(freqPeaks < auxStruct.logo_bounds(1))) % if it crosses cluster boundary was not a pure event; bounds from logothetis nature paper
                disp(tr)
      
            else  % we found a clean within cluster event 
                
                [~, mxCluster] = max(bestFactor);
                % add relevant info 
                
                % wideband power
                SNR(1:size(ripplePower, 1), :, jj)     = ripplePower;
                rippleMaxFreq(jj)                      = frequencies(mxCluster); % what frequency was this ripple best explained in
                relRipple                              = ccrf(trIndx, :); relRipple = relRipple(~isnan(relRipple)); relRipple = relRipple(relRipple~=0);
                trX                                    = relRipple(1);          
                tmpRipple                              = relRipple(4:end);
                
                if rippleMaxFreq(jj) > auxStruct.logo_bounds(3)                                    % if its not in ripple range
                         currChan(trX, (tmpRipple) + 3)         = ones(1, length(tmpRipple));
                    
                elseif (rippleMaxFreq(jj) >= auxStruct.logo_bounds(1)) && (rippleMaxFreq(jj) < auxStruct.logo_bounds(2)) % logothetis bands
                    currChanGamma(trX, (tmpRipple) + 3)         = ones(1, length(tmpRipple));
                    
                elseif rippleMaxFreq(jj) < auxStruct.logo_bounds(1) % logothetis bands
                    currChanSigma(trX, (tmpRipple) + 3)         = ones(1, length(tmpRipple)); 
                end
                
                relTrials(jj)                         = trX;
                pureRipples(jj, 1:length(relRipple))  = relRipple;
                origTrials(jj)                        = tr;
    
                jj = jj + 1;
            end
        end

        trIndx = trIndx + 1;    
    end
    
    % generate struct and save
    rippleStruct = struct;  
    
    % these we will split later on using frequencies into relevant bands
    rippleStruct.pureRipples        = pureRipples; % this is *not* only the 80hz+ ripples 
    rippleStruct.cleanTrials        = relTrials;
    rippleStruct.frequencies        = rippleMaxFreq;
    rippleStruct.SNR                = SNR;
    
    % these are already split
    rippleStruct.pureRipplesHz      = currChan;
    rippleStruct.pureRipplesGammaHz = currChanGamma;
    rippleStruct.pureRipplesSigmaHz = currChanSigma;
 
    toc
end


cleanedRipples = struct;

relReg = find(auxStruct.brain_region_channels == 5);

for i = relReg
    
    ripp = load([baseP, fl]); ripp = ripp.rippleStruct;
    nTrials = ripp.pureRipples(:, 1); nTrials = nTrials(~isnan(nTrials));

    cleanedRipples(i).pureRipplesHz = zeros(size(ripp.pureRipplesHz, 1), 2500);

    for nt = 1:length(nTrials)

        currRipp1 = ripp.pureRipples(nt, 2);
        currRipp2 = ripp.pureRipples(nt, 3);
        currRipp3 = ripp.pureRipples(nt, 4);

        if currRipp1 ~= 1
            currRipp = currRipp1;
        elseif currRipp1 == 1 & currRipp2 ~= 1
            currRipp = currRipp2;
        elseif currRipp1 == 1 & currRipp2 == 1
            currRipp = currRipp3;
        end
        cleanedRipples(i).pureRipplesHz(nTrials(nt), currRipp) = 1;
    end
end

