function [RegressorContainer] = generateAngleRegressor(AuxiliaryData)

RegressorContainer =  struct;
RegressorContainer.TrialAngle = [];
RegressorContainer.TrialVector = [];

Units = 1:length(AuxiliaryData);

angleType = 2; % 1 chosen unchosen, 2 left right; 1 doesnt have good sampling of the space

for i = Units
    tic
    
    leftCol = 1;
    rightCol = 2;
    
    % stretch it out
    unit = AuxiliaryData(i);
    strechedRank = [unit.rnkSpatVal{:}];                                      % we assume magnitude is always on the x axis of our 2d space and probability is always on the y axis
    probPos = unit.Prob_top;
    
    Trials = 1:length(probPos);
    
    chosM    = calc_rank([unit.ChML + unit.ChMR]);
    unchosM  = calc_rank([unit.UnPL + unit.UnPR]);
    
    chosP    = calc_rank([unit.ChPL + unit.ChPR]);
    unchosP  = calc_rank([unit.UnPL + unit.UnPR]);
    
    for tr = Trials
        
        if probPos(tr) == 1                                                       % 1 - top, 2 - bottom
            
            MagL =      strechedRank(2, leftCol);
            ProbL =     strechedRank(1, leftCol);
            MagR =      strechedRank(2, rightCol);
            ProbR =     strechedRank(1, rightCol);
            
        elseif probPos(tr) == 2
            
            MagL =      strechedRank(1, leftCol);
            ProbL =     strechedRank(2, leftCol);
            MagR =      strechedRank(1, rightCol);
            ProbR =     strechedRank(2, rightCol);
            
        end
        
        % left vs right
        
        if angleType == 2
            
            xLen = MagR - MagL;
            yLen = ProbR - ProbL;
            
            RegressorContainer(i).TrialVector(tr, :) = [xLen, yLen];
            
            
            Angle = atan2(yLen, xLen)*180/pi;
            
            RegressorContainer(i).TrialAngle(tr) = mod(Angle, 360);
            
            leftCol = leftCol + 2;
            rightCol = rightCol + 2;
            
        else
            % chosen unchosen
            
            xLen = chosM - unchosM;
            yLen = chosP - unchosP;
            
        end
    end
    
    if angleType == 1
        
        RegressorContainer(i).TrialVector = [xLen, yLen];
        
        
        Angle = atan2(yLen, xLen)*180/pi;
        
        RegressorContainer(i).TrialAngle = mod(Angle, 360);
        
    end
    
    % Exclude trials with xLen = 0 & yLen = 0
    toExclude = find(sum(RegressorContainer(i).TrialVector(:,:) == 0, 2) == 2);  % if the sum of a row equals to 2, it means both xLen and yLen are 0, meaning there was no diff between both on a 2d space.
    indxVec = ones(length(Trials), 1);
    indxVec(toExclude) = 0;
    indxVec = find(indxVec);
    
    
    RegressorContainer(i).TrialAngle =  RegressorContainer(i).TrialAngle;
    RegressorContainer(i).TrialVector = RegressorContainer(i).TrialVector(indxVec, :);
    RegressorContainer(i).excludeTrials = indxVec;
    
    
    [i toc];
end     % end units

end