function [RegressionResults] = run_qf_PE_Perm_glmfit(chanInfo, auxStruct, regressOut)

crosvals   = auxStruct.crosvals;
symmetries = auxStruct.symmetries;

Units = 1:length(chanInfo);
RegContainerTest = zeros(3, length(Units), length(symmetries), crosvals);
orientations = zeros(crosvals, length(Units), length(symmetries));
startBin     = round((auxStruct.preEvent + auxStruct.slSizeStart)/auxStruct.slideWidth + 1);
endBin       = round((auxStruct.preEvent + auxStruct.slSizeEnd)/auxStruct.slideWidth + 1);

for u = Units
    
    parforTmpSHALF = zeros(3, length(symmetries), crosvals);
    tmpOrientations = zeros(crosvals, length(symmetries));
    
 	FR = nanmean(chanInfo(u).Y(:, startBin:endBin), 2);
    
    trialHeading = chanInfo(u).Angle;

    trLen = size(FR, 1);
    
    % optional - regress out value and compute results on residuals 

    if regressOut
        out = regressOutValue(chanInfo, u, auxStruct, 1); % 2 SR - 1 PE
        FR  = out;
    end

    for crosval = 1:crosvals % need to make sure this is same as prev
        
        trainTr = 1:trLen;
        testTr  = 1:trLen;
                
        % this inverts the crosval 
        trainTr((crosvals + 1 - crosval):crosvals:trLen)  = [];
        testTr = setdiff(testTr, trainTr);
        
        testReg = zeros(length([trainTr testTr]), length(symmetries));
        
        trainY = zscore(FR(trainTr));
        testY  = zscore(FR(testTr));
         
        for field = 1:length(symmetries)
            
            TrainingReg = [cosd(symmetries(field)*trialHeading)', sind(symmetries(field)*trialHeading)'];
            
            [~,~,stats] = glmfit(TrainingReg(trainTr, :), trainY);
            
            betas = stats.beta(2:3);
            
            orientation = atan2d(betas(2), betas(1)) / (symmetries(field));
            
            adjustedHeading = trialHeading - orientation;
            tmpOrientations(crosval, field) = orientation;
            
            testReg(:, field) = cosd(symmetries(field)*adjustedHeading)';
            
            [~,~,stats] = glmfit(testReg(testTr, field), testY);
            
            out = [stats.beta(2), stats.t(2), stats.p(2)]; 
            parforTmpSHALF(:, field, crosval) = out';
            
        end
        
    end     % close crosval

    RegContainerTest(:, u, :, :) = parforTmpSHALF;
    orientations(:, u, :) = tmpOrientations;
    
end         % close units

RegressionResults.Test = RegContainerTest;
RegressionResults.Omegas = orientations;

end






