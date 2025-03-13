function channel_info = compute_ripple_aligned_lfp_activity()

try

    load('channel_info')

catch

    respMat   = load('resp mat lfp d');
    files = 'folder w rip fl';

    relReg   = find(auxStruct.brain_region_channels == 5);

    % we will add titles in ppt
    absIndx = 1;
    absIndxChan     = 1;

    rawRippleLFP    = nan(8000, 2000);
    spectralPower   = nan(8000, 2000, 133);
    channelAvg      = [];
    channelAvg_full = [];
    spectralPowerAvg = nan(97, 2000, 133);

    auxStruct.powerComp = 0;
    clipping_wind = 50;

    for i = relReg

        rawRippleLFP_Full = nan(8000, 2000);
        ripp = load('ripples');
        tic

        % bandpass raw signal
        currChanP                 = preprocessLFP(respMat(i).ch, 0, auxStruct.rippleBounds, auxStruct.Fs, 14);
        currChan                  = ripp.pureRipples;
        currChanLFP               = currChanP.lfp;
        currChanP                 = preprocessLFP(respMat(i).ch, 0, [4, 100], auxStruct.Fs,  14);
        currChanLFP_Full          = currChanP.lfp;                                     

        tot_len = size(currChanLFP, 2);

        % get ripples
        currChan  = currChan(~isnan(currChan(:, 1)), :);
        relFreqx  = find(ripp.frequencies > auxStruct.logo_bounds(3));
        currChan  = currChan(relFreqx, :);

        % define unique trials 
        unqTrials = unique(currChan(:, 1));

        if auxStruct.powerComp
            currPower                 = preprocessLFP(respMat(i).ch, 1, auxStruct.widebandBounds_supps, auxStruct.Fs, 20);
            currPower                 = currPower.lfpMagnitude.^2;
        end

        for tr = 1:size(unqTrials, 1)
            
            currRippleEvent               = currChan(unqTrials(tr) == currChan(:, 1), :);
            currTrialLFP                  = currChanLFP(unqTrials(tr), :);
            currTrialLFP_Full             = currChanLFP_Full(unqTrials(tr), :);

            if auxStruct.powerComp
                currTrialPower                = currPower(unqTrials(tr), :, :);
            end

            for rp = 1:size(currRippleEvent, 1)

                tmp = currRippleEvent(rp, 4:end);
                tmp = tmp(~isnan(tmp));
                tmp = tmp(tmp>0);
                completeTiming = tmp(1)-auxStruct.offset:tmp(end)+auxStruct.offset;

                if all(completeTiming > clipping_wind) && all(completeTiming < (tot_len - clipping_wind)) % avoid clipping because we're adding data on both sides

                    alignedRippleLFP = currTrialLFP(:, tmp);
                    [~,b]            = min(alignedRippleLFP);

                    dff = auxStruct.maxIndx - b - auxStruct.offset;

                    alignedRippleLFP       = currTrialLFP(:, completeTiming);
                    alignedRippleLFP_Full  = currTrialLFP_Full(:, completeTiming);
                    alignedPower           = currTrialPower(:, completeTiming, :);

                    rawRippleLFP(absIndx, (1 + dff):(length(alignedRippleLFP)+dff)) = alignedRippleLFP;
                    rawRippleLFP_Full(absIndx, (1 + dff):(length(alignedRippleLFP_Full)+dff)) = alignedRippleLFP_Full;
                 
                    if auxStruct.powerComp
                        spectralPower(absIndx, (1 + dff):(length(alignedRippleLFP_Full)+dff), :) = squeeze(alignedPower);
                    end

                    absIndx = absIndx + 1;

                else
                    continue
                end
            end
        end

        if auxStruct.powerComp
            spectralPowerAvg(absIndxChan, :, :) = squeeze(nanmean(spectralPower));
        end

        channelAvg(absIndxChan, :) = nanmean(rawRippleLFP, 1);
        channelAvg_full(absIndxChan, :) = nanmean(rawRippleLFP_Full);
        absIndxChan = absIndxChan + 1;
        toc

    end

    channel_info = struct;
    channel_info.channelAvg          = channelAvg;
    channel_info.channelAvg_full     = channelAvg_full;
    channel_info.channelAvg_spectral = spectralPowerAvg;

    % save channel_info 

end
end
