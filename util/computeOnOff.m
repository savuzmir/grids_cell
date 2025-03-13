function [out1, out2, out3] = computeOnOff(chanInfo, inputT, omegas, fold_sym, wavelength, cvals)

out2 = nan;

for i = 1:size(chanInfo, 2)
    
    Y = nanmean(chanInfo(i).Y(:, inputT), 2);    
    trialHeading = chanInfo(i).Angle;

    on_off_signal = nan(100, cvals, fold_sym*2);
    
    inv = cvals;
    
    for cv = 1:cvals
        
        part = inv:cvals:size(Y, 1);
        inv = inv - 1;
        
        currAng    = mod(trialHeading - omegas(cv, i), 360); currAng = currAng(part)';

        currY      = Y;  currY = zscore(currY(part));
        
        for peak = 1:fold_sym
            % find on and off centres:
            on = mod((peak-1)*wavelength, 360);
            off = mod((peak-1)*wavelength + wavelength/2, 360);
            % find on and off trials
            if on-wavelength/4<0
                on_trials = [find(currAng>mod(on-wavelength/4, 360)); find((currAng>on-wavelength/4).*(currAng<on+wavelength/4))];
                off_trials = find((currAng>off-wavelength/4).*(currAng<off+wavelength/4));
            else
                on_trials = find((currAng>on-wavelength/4).*(currAng<on+wavelength/4));
                off_trials = find((currAng>off-wavelength/4).*(currAng<off+wavelength/4));
            end
            on_off_signal(1:length(on_trials), cv, (peak-1)*2+1) = currY(on_trials);
            on_off_signal(1:length(off_trials), cv, (peak-1)*2+2) = currY(off_trials);
           
        end
    end
    
    all_ons = on_off_signal(:, :, 1:2:end);
    all_offs = on_off_signal(:, :, 2:2:end);
    all_ons = all_ons(:);
    all_offs = all_offs(:);
    
    [~, p, ~, stats] = ttest2(all_ons(~isnan(all_ons)), all_offs(~isnan(all_offs)));
    
    on_off_signal_log(:, :, :, i) = [on_off_signal];
    containerAligned(i, :) = [p, stats.tstat];
end

out1 = on_off_signal_log;
out3 = containerAligned;

end
