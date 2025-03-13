function residOut = regressOutValue(chanInfo, u, auxStruct, type)
% note - this function assumes both X and Y are already sorted; other vals have to be sorted 

        if type == 1
            Y = nanmean(chanInfo(u).Y(:, auxStruct.startBin:auxStruct.endBin), 2);
        else
            Y = chanInfo(u).Y; 
        end

        % first make sure alignment is correct
        x1 = chanInfo(u).Aux.Right_pay - chanInfo(u).Aux.Left_pay;
        y1 = chanInfo(u).Aux.Right_prob - chanInfo(u).Aux.Left_prob;
        angles = atan2(y1, x1)*180/pi;
        [~, tmp] = sort(mod(angles, 360));
    
        % sort data to make sure regressors are sensible 
        leftMag   = chanInfo(u).Aux.Left_pay(tmp);
        rightMag  = chanInfo(u).Aux.Right_pay(tmp);
        
        leftProb  = chanInfo(u).Aux.Left_prob(tmp);
        rightProb = chanInfo(u).Aux.Right_prob(tmp);        
        
        leftValue  = leftMag.*leftProb;
        rightValue = rightMag.*rightProb;
        
        chosenMag   = (chanInfo(u).Aux.ChML(tmp) + chanInfo(u).Aux.ChMR(tmp));
        chosenProb  = (chanInfo(u).Aux.ChPL(tmp) + chanInfo(u).Aux.ChPR(tmp));
         
        chosenValue = chosenMag.*chosenProb;
        
        % unchosen val 
        unchLeft = (chanInfo(u).Aux.ChML(tmp) > 0);
        unchRight = (chanInfo(u).Aux.ChMR(tmp) > 0);

        % what remains below is unchosen
        tmpMagR = rightMag; tmpMagL = leftMag;
        tmpProbR = rightProb; tmpProbL = leftProb;
        
        tmpMagL(find(unchLeft)) = 0; tmpProbL(find(unchLeft)) = 0;
        tmpMagR(find(unchRight)) = 0; tmpProbR(find(unchRight)) = 0;
        
        unchosenMag  = tmpMagL  + tmpMagR;
        unchosenProb = tmpProbL + tmpProbR;
        unchosenValue  = unchosenMag.*unchosenProb;
  
        % chosen side 
        chRight = unchLeft; 

        % what is on top 
        probTop = [];

        for n = 1:size(chanInfo(u).Aux.spatial_values, 2)
            probTop(n) = any(unique(chanInfo(u).Aux.Left_pay) == chanInfo(u).Aux.spatial_values{n}(1)); % if top left is equal to any of the left pay, then mag was on top
        end

        % sort correctly 
        probTop = probTop(tmp)'; 
        dm = [ones(size(probTop, 1), 1), chRight, probTop, chanInfo(u).sessRT, abs(leftValue-rightValue), chosenValue-unchosenValue, ...
             chosenMag, unchosenProb, chosenProb, unchosenMag];      
        
        [~, ~, ~, ~, resid] = ols(zscore(Y), dm, eye(size(dm, 2)));
        residOut = resid;
end