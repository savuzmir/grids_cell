function [RegressionResults] = run_qf_PE_wholeSess(chanInfo, auxStruct, trainOrientations, trainTest, regressOut)

symmetries = auxStruct.symmetries;
RegContainerTest = zeros(3, length(symmetries));
orientations = zeros(1, length(symmetries));

tmpBetas        = zeros(3, length(symmetries));
tmpOrientations = zeros(length(symmetries), 1);

% take the one that is given
FR = nanmean(chanInfo(1).Y(:, auxStruct.startBin:auxStruct.endBin), 2);
trialHeading = chanInfo(1).Angle;

trLen = size(FR, 1);

testReg = zeros(trLen, length(symmetries));

if ~trainTest % if 0 we just run the test to extract omegas from whole session
    
    if regressOut
        out = regressOutValue(chanInfo, 1, auxStruct, 1); % 2 SR - 1 PE
        FR = out;
    end

    % we z score the signal    
    trainY = zscore(FR);


for field = 1:length(symmetries)
    
    TrainingReg = [cosd(symmetries(field)*trialHeading)', sind(symmetries(field)*trialHeading)'];
    [~,~,stats] = glmfit(TrainingReg, trainY);
    betas = stats.beta(2:3);
    orientation = atan2d(betas(2), betas(1)) / (symmetries(field));
    tmpOrientations(field) = orientation;
    
end % close field

else
    
    if regressOut
        out = regressOutValue(chanInfo, 1, auxStruct, 1); % 2 SR - 1 PE
        FR = out;
    end

    testY  = zscore(FR);
    
    for field = 1:length(symmetries)
        
        adjustedHeading = trialHeading - trainOrientations(field);
        testReg(:, field) = cosd(symmetries(field)*adjustedHeading);
        [~, ~, stats] = glmfit(testReg(:, field), testY);
        tmp = [stats.beta(2); stats.t(2); stats.p(2)];
        tmpBetas(:, field) = tmp;
        
    end
    
end         % close units

RegContainerTest(:, :) = tmpBetas;
orientations(1, :) = tmpOrientations';

if trainTest % if we tested we want the betas
    RegressionResults.Test = RegContainerTest;
else         % if we trained we want the ori
    RegressionResults.Omegas = orientations;
end

end






