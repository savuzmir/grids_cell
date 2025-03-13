function [out1, testBetasCells, channelOrientationsFull, trainOrientations] = run_LFP_to_UNIT_QF_PE(chanInfo, unitInfo, auxStruct, regressOut)

nPerm        = 1; % we dont want to run permutations here auxStruct.nPerm;
symmetries   = auxStruct.symmetries;
unitSessions = auxStruct.unitSessions;

trainTest = 0;
%--------------------------------------
% train on LFP
%--------------------------------------
trainOrientations= [];
nChans      = size(chanInfo, 2);
nUnits      = size(unitInfo, 2);

for i = 1:nChans
    [RegressionResults]       = run_qf_PE_wholeSess(chanInfo(i), auxStruct, [], trainTest, regressOut);
    trainOrientations(i, :)   = RegressionResults.Omegas;
end

disp('trained')

% we upsample the channel betas to cell betas
channelOrientationsFull  = [];

for i = 1:nChans
    
    currChan = find(round(unitSessions(:, 4)) == i);
    repOrientations = repmat(trainOrientations(i, :), [size(currChan, 1), 1]);
    
    if currChan
        channelOrientationsFull(currChan, :) = repOrientations;
    else
        continue
    end
end

disp('upsampled')

%--------------------------------------
% test on cells
%--------------------------------------
regressOut = 0; % already regressed out so should be 0 here

tmp = struct;
try
    trainTest = 1;

    testBetasCells = nan(nUnits, length(symmetries), nPerm);

    for np = 1:nPerm
        tic
        if np > 1

            for u = 1:nUnits
                tmp(u).Y        = unitInfo(u).Y(shuffle(1:size(unitInfo(u).Y, 1)), :);
                tmp(u).Angle    = unitInfo(u).Angle;
            end

        else
            for u = 1:nUnits
                tmp(u).Y        = unitInfo(u).Y;
                tmp(u).Angle    = unitInfo(u).Angle;
            end
        end

        for j = 1:size(unitInfo, 2) 

            [RegressionResults] = run_qf_PE_wholeSess(tmp(j), auxStruct, channelOrientationsFull(j, :), trainTest, regressOut);
            testBetasCells(j, :, np)   = RegressionResults.Test(1, :);

        end
        toc
    end
    disp('tested')
catch
    ...
end
out1 = nan;
end