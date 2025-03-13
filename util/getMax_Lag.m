function [mainOut, toExcludeTrials] = getMax_Lag(mainIn, auxStruct, field, summary)
% this function takes mainIn (STRUCT) with raw or binned lfp data (nTrials x Time) in field (STRING) and returns the lag + max peak of the data
% looking at a window defined in auxStruct (STRUCT) from startSearchAt (INT) to endSearchAt (INT)
% if summary (INT) is set to 1, returns trials to exclude due to noise 

Units = 1:size(mainIn, 2);

eventOnset = auxStruct.onset;                                                  % defines search window within actual raw matrix
startSearchAt = eventOnset + auxStruct.startSearch;                            
endSearchAt = auxStruct.endSearch;                                             

toExcludeTrials = struct;

cutOff = length(startSearchAt:endSearchAt);
cutOffCrit = 0.1;  % was 0.1                                                              
cutOff = round(cutOff*cutOffCrit);

if ~summary
    mainOut = struct;
    mainOut(1).rawValMax = [];
    mainOut(1).lagMax = [];
    mainOut(1).rawValMin = [];
    mainOut(1).lagMin = [];
    
    for u = Units
        
        trials = 1:size(mainIn(u).(field)(:, startSearchAt:endSearchAt), 1);
        for tr = trials

            unitSlice = mainIn(u).(field)(tr, startSearchAt:endSearchAt);
            [rawValMax, lagMax] = max(unitSlice);                                
            [rawValMin, lagMin] = min(unitSlice);
            
            lagMax = lagMax + eventOnset; lagMin = lagMin + eventOnset;           
            
            mainOut(u).rawValMax(tr)    = rawValMax;
            mainOut(u).lagMax(tr)       = lagMax;
            mainOut(u).rawValMin(tr)    = rawValMin;
            mainOut(u).lagMin(tr)       = lagMin;
            
        end
    end
    
else
    
    for u = Units
        
        unitSlice = mainIn(u).(field)(:, startSearchAt:endSearchAt);
        avgUnitSlice = nanmean(unitSlice);
%       avgUnitSlice = mean(unitSlice);
        [rawValMax, lagMax] = max(avgUnitSlice);                                 
        [rawValMin, lagMin] = min(avgUnitSlice);
        

        lagMax = lagMax + eventOnset; lagMin = lagMin + eventOnset;              
        
        mainOut(u, :) = [rawValMax, lagMax, rawValMin, lagMin];                    

        ts_mean = nanmean(unitSlice);
        ts_std = nanstd(unitSlice);
        
%       ts_mean = mean(unitSlice);
%       ts_std  = std(unitSlice);

        ub = ts_mean + 3*ts_std;
        lb = ts_mean - 3*ts_std;
        
        r = size(unitSlice', 2);
        
        ub = repmat(ub', [1, r]);
        lb = repmat(lb', [1, r]);
        
        excludeTrsUB = unitSlice' > ub;
        excludeTrsLB = unitSlice' < lb;
        
        alltr = excludeTrsUB + excludeTrsLB; % cant be too high and too low at the same time 
        trials = find(sum(alltr) > cutOff);
        
        nTrExcl = length(trials)/size(unitSlice, 1)*100;
        
        toExcludeTrials(u).trialNums = trials;
        toExcludeTrials(u).trialsPerc = nTrExcl;

        
    end
    
end

end






