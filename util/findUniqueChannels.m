function [correctCells, sessions, all_info, origIndx] = findUniqueChannels(allFiles, all_units)
% Compute the actual channels upon which we are going to do the regressions over - done for the sim trials
% this requires a structure with subject, session number and channel number information

Units = 1:length(allFiles);

prevSess = 1;

chanID = 1;
all_info = [];

for u = Units
    
    currSub =  all_units{u}.UnitInfo.subject;
    currSess = all_units{u}.UnitInfo.session_number;
    currChan = all_units{u}.UnitInfo.channel_number;
    currAnimal = [all_units{u}.UnitInfo.isfrank, all_units{u}.UnitInfo.ismiles]; 
    
    if u >= 2
        if prevSess ~= currSess                                              % needs to be separate to evaluate properly
            chanID = chanID + 1;
        else
            if prevChan ~= currChan
                chanID = chanID + 1;
            end
        end
        
    end
    
    prevSess = currSess;
    prevChan = currChan;
    
    all_info(u,:) = [u, currSess, currChan, chanID, currAnimal];
end

[~, origIndx] = unique(all_info(:,4));

 correctCells = all_info(origIndx, 1);                                  % these are the channels we are going to index now
 sessions = all_info(origIndx, 2);
end
