function figure1_BehaviourValue()
clear; clf;

%% load data structs
% cd(go into data folder) 

rel_files = {'auxStruct', 'auxLFP', 'AuxiliaryCells', 'allUnits'};
for fl = 1:length(rel_files)
    load(rel_files{fl});
end

%% compute accuracy for each chan
accuracyOverall = [];

for i = 1:auxStruct.num_channels
    choseLeft = auxLFP(i).ch.ChML > 0;
    EV = (auxLFP(i).ch.Left_pay.*auxLFP(i).ch.Left_prob) > (auxLFP(i).ch.Right_pay.*auxLFP(i).ch.Right_prob);
    accuracyOverall(i)  = sum(EV==choseLeft)/length(EV); 
end

% reduce to unique sessions
indx = 1;
avgAcc = [];

for an = 0:1
    
    currSessions = auxStruct.chanSessions(auxStruct.chanSessions(:, 6) == an, 2);
    unqSessions = unique(currSessions);

    for us = unqSessions'
        avgAcc(indx) = squeeze(nanmean(accuracyOverall(currSessions == us)));
        indx = indx + 1;
    end
end

%% compute psychometric curve

leftMagAll      = []; 
leftProbAll     = [];
rightMagAll     = [];
rightProbAll    = [];
chosenAll       = [];

relSessions = [];
sessIdx     = 1;

for i = 1:auxStruct.num_channels

    % discretize
    leftMag  = discretize(auxLFP(i).ch.Left_pay, 5);
    rightMag = discretize(auxLFP(i).ch.Right_pay, 5);
    
    leftProb  = discretize(auxLFP(i).ch.Left_prob, 5);
    rightProb = discretize(auxLFP(i).ch.Right_prob, 5);
    
    chosen    = auxLFP(i).ch.ChML == 0;
   
    % reduce tot sessions
    if i == 1 % if first of what we are interested in
        sess = auxStruct.chanSessions(i, 2);
       
        leftMagAll   = [leftMagAll; leftMag];
        leftProbAll  = [leftProbAll; leftProb];
        rightMagAll  = [rightMagAll; rightMag];
        rightProbAll = [rightProbAll; rightProb];
        chosenAll    = [chosenAll; chosen];
        
    elseif sess == auxStruct.chanSessions(i, 2)
        sess = auxStruct.chanSessions(i, 2);
        continue
    else
        sess = auxStruct.chanSessions(i, 2);

        leftMagAll   = [leftMagAll; leftMag];
        leftProbAll  = [leftProbAll; leftProb];
        rightMagAll  = [rightMagAll; rightMag];
        rightProbAll = [rightProbAll; rightProb];
        chosenAll    = [chosenAll; chosen];
        
    end
    relSessions(sessIdx) = sess;
    sessIdx = sessIdx + 1;
end

% compute curve 
left_mag   = leftMagAll;
left_prob  = leftProbAll;
right_mag  = rightMagAll;
right_prob = rightProbAll;
chosen     = chosenAll;

[out1, ~, out3] = computePsychometricCurve(left_mag, left_prob, right_mag, right_prob, chosen);

ffit1 = fitPsyche(out3, out1, 'GLM');

% get coeffs
a = ffit1.coeffs(1); b = ffit1.coeffs(2);

%% compute value 
% prepare raw rasters

% allocate info
tp_len  = length(auxStruct.timepoints_sr);
reg_len = 2;

cpdStats      = nan(auxStruct.num_cells, tp_len, reg_len, auxStruct.nPerm + 1);
cpdStatsPE    = nan(auxStruct.num_cells, reg_len, auxStruct.nPerm + 1);

% run analysis for chosen value + chosen value difference
parfor nPerm = 1:(auxStruct.nPerm + 1)
    tic
    tmpStats   = nan(auxStruct.num_cells, tp_len, reg_len);
    tmpStatsPE = nan(auxStruct.num_cells, reg_len);

    for unit = 1:size(allUnits, 2)
        
        if nPerm == 1
            Y               = allUnits(unit).trials;
            Y_PE            = allUnits(unit).trialsPE;
        else
            Y               = allUnits(unit).trials;
            Y               = Y(randperm(size(Y, 1)), :);
            Y_PE            = allUnits(unit).trialsPE(randperm(size(Y, 1)));
        end
        
        chosenVal       = [AuxiliaryCells(unit).ChEVL + AuxiliaryCells(unit).ChEVR]; 
        unchosenVal     = [AuxiliaryCells(unit).UnEVL + AuxiliaryCells(unit).UnEVR];
        chDiff          = chosenVal - unchosenVal;

        dm_ch_val = [ones(length(chosenVal), 1), chosenVal];
        [beta] = cpd(Y, dm_ch_val);
        tmpStats(unit, :, 1) = beta(2, 1:tp_len); 
        
        % compute over time
        % compute ch val 
        dm_ch_val_diff = [ones(length(unchosenVal), 1), chDiff];

        % compute ch val diff
        [beta] = cpd(Y, dm_ch_val_diff);
        tmpStats(unit, :, 2) = beta(2, 1:tp_len);
  
        % compute on PE
        [beta] = cpd(Y_PE, dm_ch_val);
        tmpStatsPE(unit, 1) = beta(2);
        [beta] = cpd(Y_PE, dm_ch_val_diff);
        tmpStatsPE(unit, 2) = beta(2);  
            
    end

    [nPerm, toc]
    cpdStats(:, :, :, nPerm)      = tmpStats;
    cpdStatsPE(:, :, nPerm)       = tmpStatsPE;
end

% generate null
num_br_reg  = length(auxStruct.areaCodes);
nullStatsPE = nan(auxStruct.nPerm, reg_len, num_br_reg);
trueStatsPE = nan(reg_len, num_br_reg);

for rg = 1:2 % 1 - chosen val, 2 - chosen val diff
    for np = 2:(auxStruct.nPerm + 1)
        for br = 1:length(auxStruct.areaCodes) % brain regs
            brag = find(auxStruct.brain_region_cells == auxStruct.areaCodes(br));
            nullStatsPE(np-1, rg, br)  = nanmean(cpdStatsPE(brag, rg, np), 1);            
        end
    end
    
    for br = 1:4
        brag = find(auxStruct.brain_region_cells == auxStruct.areaCodes(br));
        trueStatsPE(rg, br)  =  nanmean(cpdStatsPE(brag, rg, 1));
    end
end

nullStatsMean = squeeze(nanmean(nullStatsPE));
nullStatsSTD  = squeeze(nanstd(nullStatsPE));

zscorePE = (trueStatsPE - nullStatsMean)./nullStatsSTD;
normcdf(-zscorePE) % pvalues

cpdValues = nan(200, reg_len, num_br_reg);

for br = 1:length(auxStruct.relRegs)
    currReg = find(auxStruct.brain_region_cells == auxStruct.areaCodes(br));

    % ch val
    currCPD = cpdStatsPE(currReg, 1, 1);
    cpdValues(1:length(currReg), 1, br) = currCPD;
    
    % ch diff
	currCPD = cpdStatsPE(currReg, 2, 1);
    cpdValues(1:length(currReg), 2, br) = currCPD;
    
end

% start figure 
curr = figure(1);

threshUPPE = squeeze(prctile(nullStatsPE, 95, 1));

subplot(2, 2, 1)
raincloud_plot(avgAcc*100, 'box_on', 1, 'color', auxStruct.greenVec, 'bandwidth', 1, 'summary_type', 2, 'box_dodge', 1, 'cloud_edge_col', auxStruct.greenVec);
figElements(curr, '', 'Accuracy [%]', '', [50 100], [], [], {}, [], {}, 20, [], 1); box off
view([270, 90])

subplot(2, 2, 2)
figInfo = {'', 'Rank Difference (Right-Left)', 'P(choose right)', [-24, 24], [], [-24:6:24], {'-4', '-3', '-2', '-1', '0', '1', '2', '3', '4'}, [], {}, 20, [], []};       
plotPsychometricFunction(ffit1, curr, figInfo, auxStruct.greenVec); hold on

% ch ev
rel_reg = 1;

subplot(2, 2, 3)
for i = 1:length(auxStruct.relRegs)
    currReg = find(auxStruct.brain_region_cells == auxStruct.areaCodes(i));
    
    h(i) = plotmse(squeeze(cpdStats(currReg, :, rel_reg))*100, auxStruct.areaCol(i, :), [0, 1]); hold on
end
figElements(curr, 'Chosen Value Signal', 'Time [msec]', 'CPD [%]', [1 141], [0, 2.5], [24, 41 81 121], ...
                {'Cue ON', '-400', 'Choice ON', '400'}, [], {}, 20, [], []); box off

% add line for mean RT 
meanRT = [];
for i = 1:auxStruct.num_channels
    meanRT(i) = nanmean(auxStruct.sessRT{i});
end

meanRT_tot = round(nanmean(meanRT)); % 556 corresponds to ~24
plot([24, 24], [0, .05], 'color', 'k', 'linewidth', 3.5)
legend(h(:), 'ACC', 'DLPFC', 'OFC', 'VMPFC', 'Location', 'northwest', 'NumColumns', 2); legend box off

subplot(2, 2, 4)
b = bar(nanmean(squeeze(cpdValues(:, rel_reg, :)))*100, 'edgeColor', 'none'); hold on; box off 
b.FaceColor = 'flat';
for j = 1:length(auxStruct.areaCol)
    b.CData(j, :) = auxStruct.areaCol(j, :);
end

x_n = [0.5 1.5; 1.5, 2.5; 2.5, 3.5; 3.5, 4.5];
y_n = [];

for breg = 1:4
    y_n(breg, :) = [threshUPPE(rel_reg, breg)*100, threshUPPE(rel_reg, breg)*100];
end

for i = 1:4
    plot(x_n(i, :), y_n(i, :), 'k', 'linewidth', 2); hold on
end

err = nanstd(squeeze(cpdValues(:, rel_reg, :)), [], 1)./sqrt(sum(~isnan(squeeze(cpdValues(:, rel_reg, :))), 1));
errorbar(1:4, nanmean(squeeze(cpdValues(:, rel_reg, :))) * 100, err * 100,  '.', 'CapSize', 0, 'linewidth', 2, 'color', [155/255, 155/255, 155/255])
figElements(curr, '', 'Brain Regions', 'CPD [%]', [], [0, 3], [1 2 3 4], {'ACC', 'DLPFC', 'OFC', 'VMPFC'}, [], {}, 20, [], [])

%% supplementary fig 
% ch val diff

curr = figure(6)
rel_reg = 2;

subplot(1, 2, 1)
for i = 1:length(auxStruct.relRegs)
    currReg = find(auxStruct.brain_region_cells == auxStruct.areaCodes(i));
    
    h(i) = plotmse(squeeze(cpdStats(currReg, :, rel_reg))*100, auxStruct.areaCol(i, :), [0, 1]); hold on
end
figElements(curr, 'Chosen Value Difference Signal', 'Time [msec]', 'CPD [%]', [1 141], [0, 2.5], [41 81 121], ...
                {'-400', 'Choice ON', '400'}, [], {}, 20, [], []); box off
            
legend(h(:), 'ACC', 'DLPFC', 'OFC', 'VMPFC', 'Location', 'northwest', 'NumColumns', 2); legend box off

subplot(1, 2, 2)
b = bar(nanmean(squeeze(cpdValues(:, rel_reg, :))) * 100, 'edgeColor', 'none'); hold on; box off 
b.FaceColor = 'flat';
for j = 1:length(auxStruct.areaCol)
    b.CData(j, :) = auxStruct.areaCol(j, :);
end

x_n = [0.5 1.5; 1.5, 2.5; 2.5, 3.5; 3.5, 4.5];
y_n = [];

for breg = 1:4
    y_n(breg, :) = [threshUPPE(rel_reg, breg)*100, threshUPPE(rel_reg, breg)*100];
end
for i = 1:4
    plot(x_n(i, :), y_n(i, :), 'k', 'linewidth', 2); hold on
end

err = nanstd(squeeze(cpdValues(:, rel_reg, :)), [], 1)./sqrt(sum(~isnan(squeeze(cpdValues(:, rel_reg, :))), 1));
errorbar(1:4, nanmean(squeeze(cpdValues(:, rel_reg, :))) * 100, err * 100,  '.', 'CapSize', 0, 'linewidth', 2, 'color', [155/255, 155/255, 155/255])
figElements(curr, '', 'Brain Regions', 'CPD [%]', [], [0, 2.5], [1 2 3 4], {'ACC', 'DLPFC', 'OFC', 'VMPFC'}, [], {}, 20, [], [])
