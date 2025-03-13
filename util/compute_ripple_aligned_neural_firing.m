function rasterMatrixIndiv = compute_ripple_aligned_neural_firing(auxStruct)

try 
       
    % load rasterMatrixIndiv

catch
    respMat   = load('resp lfp data'); % note due to task structure, the response locked mat effectively captures most of fix / rew period 
    respMatFR = load('resp fr mat'); 
    fixMatFR  = load('fix fr mat'); 
    
    % get ripple folder info 
    rippleFolder = '';
    
    % generate baseline firing rate from first 300 msec of fixation before cue on
    baselineFR = {};
    for i = 1:auxStruct.num_cells
        baselineFR{i} = nanmean(fixMatFR(i).FixMatrix(:, auxStruct.timeSlice), 2); % average
    end
    
    % pregenerate containers for vmpfc units
    rasterMatrixIndiv = struct;
    rasterMatrixIndiv.trial = [];
    
    % preallocate lengths to avoid zeros
    container_max_num = 1800;
    container_max_len = 350;

    unitsVMPFC = length(find(auxStruct.brain_region_cells == 5));
    
    for i = 1:unitsVMPFC
        rasterMatrixIndiv(i).container = nan(container_max_len, container_max_num);
        rasterMatrixIndiv(i).container_NR = nan(container_max_len, container_max_num);
        rasterMatrixIndiv(i).baselineFR = nan(container_max_len, 1); 
    end
    
    % fetch channels
    baseP = rippleFolder;
    files = dir(baseP);
    files = {files.name};
    files = files(3:end);
    
    % rearrange to match correctly
    idx = [];
    for i = 1:length(files)
        tmp = files{i}; tmp = split(tmp, '_'); tmp = split(tmp{2}, '.');
        idx(i) = str2num(tmp{1});
    end
    
    [~, xd] = sort(idx);
    files = files(xd);
    
    chanIndx = 1;
    cellIndx = 1;
   
    containerIndx1=1;

    % run ripple locking
    for chan = find(auxStruct.brain_region_channels == 5)
        tic
    
        % filter data;
        respLFP   = preprocessLFP(respMat(chan).ch, 0, [auxStruct.rippleBounds(1), auxStruct.rippleBounds(2)], auxStruct.Fs, 14);
        tot_len = size(RespLFP.lfp, 2);
%        ripp = load('indiv ripple files');
        
        % prepare ripple events
        currChanRipple                  = ripp.pureRipples; currChanRipple = currChanRipple(~isnan(currChanRipple(:, 1)), :);
        currChanRipple = currChanRipple(~isnan(currChanRipple(:, 1)), :);
        relFreqx  = find(ripp.frequencies > auxStruct.rippleBounds(1));
        currChanRipple  = currChanRipple(relFreqx, :);
    
        % check how many trials there were
        ripple_trials = currChanRipple(:, 1);
    
        % due to how data is stored, we need to adjust the columns to gete correct information on a trial-by-trial basis
        currChanRippleSifted = nan(300, 350);
    
        for tr = 1:size(currChanRipple, 1)
        
            currRipp1 = currChanRipple(tr, 2);
            currRipp2 = currChanRipple(tr, 3);
        
            % need to fix depending on how its saved
            if currRipp1 ~= 1
                currChanRippleSifted(tr, 1:size(currChanRipple(tr, 3:end), 2)) = currChanRipple(tr, 3:end);   % we dont take one immediately after but the one after that to account for the peak freq 
            elseif currRipp1 == 1 & currRipp2 ~= 1
                currChanRippleSifted(tr, 1:size(currChanRipple(tr, 4:end), 2)) = currChanRipple(tr, 4:end);    
            elseif currRipp1 == 1 & currRipp2 == 1
                currChanRippleSifted(tr, 1:size(currChanRipple(tr, 5:end), 2)) = currChanRipple(tr, 5:end);     
            end
    
        end
    
        currChanRipple = currChanRippleSifted;
    
        % these are the units from the corresponding channel
        unitIDs = auxStruct.unitSessions(auxStruct.unitSessions(:, 4) == chan, 1);
    
        for u = unitIDs'
    
            trIndx = 1;
            currUnit = baselineFR{u};
    
            allTrials = 1:size(respMatFR(u).ResponseMatrix, 1);
    
            no_ripple_trials = setdiff(allTrials, ripple_trials);
            controlContainerIdx = 1;
    
            for tr = no_ripple_trials          % generate control fr with no ripples during period when highest fr expected in the null 
    
                for peak = 1:size(auxStruct.peaks, 1)
    
                    % take the ripple band info
                    lfpCurr = respLFP.lfp(tr, auxStruct.peaks(peak, :));
    
                    % find the min peak of the oscillation
                    [~, indx] = min(lfpCurr);
    
                    % we define windows around and create temporal container to store fr 
                    rpLen  = length(auxStruct.peaks(peak, :)); maxWindow = 600 - rpLen; 
                    lb = floor(maxWindow/2); ub = ceil(maxWindow/2);         
                    trCurr = (auxStruct.peaks(peak, 1) - lb):(auxStruct.peaks(peak, end) + ub);  % ensures equal length
                    modif = auxStruct.maxIndx - indx - lb; % align to indx
    
                    rasterMatrixIndiv(cellIndx).container_NR(controlContainerIdx, (1+modif):(modif + size(trCurr, 2))) = respMatFR(u).ResponseMatrix(tr, trCurr);
                    controlContainerIdx = controlContainerIdx + 1;
                    
                end
            end
    
            % for trials where there was a recorded ripple
            for tr = ripple_trials'
                
                    % take current trial and from the current trial remove all nans
                    trCurr = currChanRipple(trIndx, :);
                    trCurr = trCurr(~isnan(trCurr)); 
    
                    % remove zeroes in case these were added 
                    trCurrOrig = trCurr; trCurrOrig = trCurrOrig(trCurrOrig > 0);
                    trCurr = trCurr(trCurr > 0);
    
                    % we generate a window around the ripple such that we
                    rpLen  = length(trCurr); maxWindow = 350 - rpLen;
                    lb = floor(maxWindow/2); ub = ceil(maxWindow/2); % we make it such that we add at the end
                    trCurr = (trCurrOrig(1) - lb):(trCurrOrig(end) + ub); % this will make all of them have the same length
                    
                    if (any(trCurr < 1) | any(trCurr > tot_len)) % prevent clipping
                        ... 
                    else
    
                        % take the ripple band info
                        lfpCurr = respLFP.lfp(tr, trCurrOrig);
    
                        % find the lowest dip of the oscillation
                        [~, indx] = min(lfpCurr);
    
                        modif = auxStruct.maxIndx - (indx + lb); 
                        rasterMatrixIndiv(cellIndx).baselineFR(trIndx, :) = currUnit(tr);
            
                        rasterMatrixIndiv(cellIndx).container(containerIndx1, (1+modif):(modif + size(trCurr, 2))) = respMatFR(u).ResponseMatrix(tr, trCurr);
                        rasterMatrixIndiv(cellIndx).baselineFR(containerIndx1) = currUnit(tr);
                        
                        containerIndx1 = containerIndx1 + 1;
                    end
                
                trIndx = trIndx + 1; % this adds nans so need to remove later 
            end
    
            cellIndx = cellIndx + 1;
            containerIndx1 = 1;
    
        end
        toc
        chanIndx = chanIndx + 1;
    
    end
    
    % save rasterMatrixIndiv

end