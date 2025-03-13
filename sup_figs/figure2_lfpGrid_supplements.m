function figure2_lfpGrid_supplements(auxStruct)
clear; clf;

%% load data structs
% cd(go into data folder) 

rel_files = {'auxStruct', 'auxLFP', 'out_results_d1', 'freq_supplements', 'out_results_d2'};
for fl = 1:length(rel_files)
    load(rel_files{fl});
end

%%
n_regs = length(auxStruct.areaCodes);

for rg = 1:n_regs
    currRg = auxStruct.relRegs{rg};

    % plot averaged stats
    avgBetasPE  = [];
    avgBetasSR  = [];
    avgSignalPE  = [];
    dy = [];

    indx = 1;
    sessions = auxStruct.chanSessions(currRg, :);

    for an = 0:1

        currSessions = sessions(sessions(:, 6) == an, 2);
        tt = auxStruct.dayInformation_stimSet(sessions(:, 1));
        tt = tt(sessions(:, 6) == an);

        unqSessions = unique(currSessions);
        rlReg = intersect(currRg, sessions(sessions(:, 6) == an, 4));

        betasPE = squeeze(nanmean(out_dat1.RegressionResults_glmfit.Test(1, rlReg, :, :), 4));
        betasSR = squeeze(nanmean(out_dat1.TestReal(rlReg, :, :, :), 3));
        avgSignal = out_dat1.signal(:, rlReg);


        for us = unqSessions'
            avgBetasPE(indx, :)    = nanmean(betasPE(currSessions == us, :), 1);
            avgBetasSR(indx, :, :) = squeeze(nanmean(betasSR(currSessions == us, :, :), 1));
            avgSignalPE(:, indx)    = squeeze(nanmean(avgSignal(:, currSessions == us), 2));
            dy(indx)               = median(tt(currSessions == us));
            indx = indx + 1;
        end
    end

    sesBetasFullPE{rg} = avgBetasPE;
    sesBetas{rg}       = avgBetasPE(:, 3);
    sesBetasSR{rg}     = squeeze(avgBetasSR(:, 3, :));
    sesBetasFullSR{rg} = avgBetasSR;
end


%% supplements 

confInt = [0, 0];
sixfold = nan(41, n_regs); % 41 is max nr sess for prealloc

for i = 1:n_regs
    sixfold(1:length(sesBetas{i}), i) = sesBetas{i}';
end

curr=figure(116)

% comparison of brain regions 
subplot(5, 3, 1)
figInfo = {'', '', 'Betas [a.u.]', [], [], [1 2 3 4], ...
    {'ACC', 'DLPFC', 'OFC', 'VMPFC'}, [], {}, 24, [], []};
[curr] = plotBar(sixfold, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, confInt, 0);

[~,b,~,d]= ttest2(sixfold(:, 4), sixfold(:, 1))
% 0.0893
[~,b,~,d]= ttest2(sixfold(:, 4), sixfold(:, 2))
% 0.0374
[~,b,~,d]= ttest2(sixfold(:, 4), sixfold(:, 3))
% 0.1077
[~,b,~,d]= ttest2(sixfold(:, 3), sixfold(:, 1))
% 0.8725

% all vmpfc periodicities over time
subplot(5, 3, 2) 

onsetModif = 0;
onsetBin = 81 - onsetModif;
endModif   = 141;
nTimepoints = length(1:(endModif-onsetModif));

regBetas = {};

rectangle('Position', [51, -0.1, 30, 0.215], ...
    'FaceColor', [1/255, 114/255, 178/255, 0.2], ...
    'EdgeColor', 'none'); hold on

plot(1:(endModif-onsetModif), repmat(0, [1, nTimepoints]), '--k', 'linewidth', 1.5); hold on
plot(repmat(onsetBin, 2, 1), [-0.1, -0.09], 'linewidth', 5, 'color', 'k');

relBetasSR = sesBetasFullSR{4};

for i = [1 2 4 5 3]
    if i ~= 3
        plotmse(squeeze(relBetasSR(:, i, :)), auxStruct.grayVec, [0, 1], 0.5)
    else
        h = plotmse(squeeze(relBetasSR(:, i, :)), auxStruct.greenVec, [0, 1], 0.6)
    end
end

figElements(curr, '', 'Time [msec]', 'Betas [a.u.]', [1 141], [-0.1, 0.115], [41-onsetModif 81-onsetModif 121-onsetModif], ...
                {'-400', 'Choice ON', '400'}, [], {}, 24, [], []); box off

% variance
subplot(5, 9, 37:39) 
TestOmegas = out_dat1.TestOmegas;
omegas = squeeze(permute(squeeze(deg2rad(TestOmegas(:, :, :, :))), [1 3 4 2]));

symmetries = 1:length(auxStruct.symmetries);

avgedOmegas = [];
stdOmegas   = [];

for u = 1:auxStruct.num_channels
    for field = symmetries
        for t = 1:121
            currOmegas = omegas(u, field, t, :);
            
            fMed = circ_std(squeeze(currOmegas));
            fVar = circ_std(squeeze(currOmegas));
            
            avgedOmegas(u, t, field) = fMed;
            stdOmegas(u, t, field) = fVar;

        end
    end
end

for rg = 1:n_regs
    currRg = auxStruct.relRegs{rg};
    
    % plot averaged stats
    avgBetasSR  = [];
    
    indx = 1;
    sessions = auxStruct.chanSessions(currRg, :);
    
    for an = 0:1
       
        currSessions = sessions(sessions(:, 6) == an, 2);
        unqSessions = unique(currSessions);
        rlReg = intersect(currRg, sessions(sessions(:, 6) == an, 4));
       
        betasSR = squeeze(avgedOmegas(rlReg, :, :));
        for us = unqSessions'
            avgBetasSR(indx, :, :) = squeeze(circ_mean(betasSR(currSessions == us, :, :))); 
            indx = indx + 1;
        end 
    end
    
    sesBetasSR{rg} = squeeze(avgBetasSR(:, :, :));
end

for rg = 1:n_regs
    if rg == 1
       rectangle('Position', [51, 4, 30, 10], ...
       'FaceColor', [1/255, 114/255, 178/255, 0.2], ...
       'EdgeColor', 'none'); hold on
       
    end
    
    if rg == 1
        plot(repmat(81, 2, 1), [4, 4.25], 'linewidth', 5, 'color', 'k');
    end
    
    mn = circ_mean(squeeze(sesBetasSR{rg}));
    lb_se = mn - circ_std(squeeze(sesBetasSR{rg}))./sqrt(size(squeeze(sesBetasSR{rg}), 1));
    ub_se = mn + circ_std(squeeze(sesBetasSR{rg}))./sqrt(size(squeeze(sesBetasSR{rg}), 1));
    
    
    mn = squeeze(mn(:, :, 3));
    lb_se = squeeze(lb_se(:, :, 3));
    ub_se = squeeze(ub_se(:, :, 3));
    
    plotmseSumm((1:121)', rad2deg([mn', lb_se', ub_se']), auxStruct.areaCol(rg, :));
    
end

figElements(curr, '', 'Time [msec]', 'Grid orientation SD [deg]', [1 121], [], [41 81 121], ...
    {'-400', 'Choice ON', '400'}, [], {}, 24, [], []); box off

% get sixfold
acc    = sesBetasSR{1}; acc = acc(:, :, 3);
dlpfc  = sesBetasSR{2}; dlpfc = dlpfc(:, :, 3);
ofc    = sesBetasSR{3}; ofc = ofc(:, :, 3);
vmpfc  = sesBetasSR{4}; vmpfc = vmpfc(:, :, 3);

[a,b,c,d]=ttest2(nanmean(acc(:, auxStruct.PE_roi), 2), nanmean(vmpfc(:, auxStruct.PE_roi), 2))
[a,b,c,d]=ttest2(nanmean(dlpfc(:, auxStruct.PE_roi), 2), nanmean(vmpfc(:, auxStruct.PE_roi), 2))
[a,b,c,d]=ttest2(nanmean(ofc(:, auxStruct.PE_roi), 2), nanmean(vmpfc(:, auxStruct.PE_roi), 2))

% different frequencies
[~,~,~, stats] = ttest(sesBetasFREQ);

subplot(5, 3, 3) 
figInfo = {'', 'Frequencies', 'Betas [a.u.]', [], [-0.02, 0.13], [1 2 3 4], {'Theta', 'Beta', 'Low Gamma', 'High Gamma'}, [], {}, 24, [], []};
[curr] = plotBar(sesBetasFREQ, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, confInt, 0);

% indiv regions 
subplot(5, 3, 4)
figInfo = {'ACC: sessions', 'Symmetries', 'Betas [a.u.]', [], [-0.06, 0.12], [1 2 3 4 5], {'4', '5', '6', '7', '8'}, [], {}, 24, [], []};
[curr] = plotBar(sesBetasFullPE{1}, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, confInt, 0);
subplot(5, 3, 5)
figInfo = {'DLPFC: sessions', 'Symmetries', 'Betas [a.u.]', [], [-0.06, 0.12], [1 2 3 4 5], {'4', '5', '6', '7', '8'}, [], {}, 24, [], []};
[curr] = plotBar(sesBetasFullPE{2}, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, confInt, 0);
subplot(5, 3, 6)
figInfo = {'OFC: sessions', 'Symmetries', 'Betas [a.u.]', [], [-0.06, 0.12], [1 2 3 4 5], {'4', '5', '6', '7', '8'}, [], {}, 24, [], []};
[curr] = plotBar(sesBetasFullPE{3}, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, confInt, 0);

% across days
betas_PE = sesBetasFullPE{4};

allBetas1   = betas_PE(dy==1 ,:);
allBetas234 = betas_PE(dy==2 | dy==3 | dy==4 ,:);

allSignal1   = avgSignalPE(:, dy==1);
allSignal234 = avgSignalPE(:, dy==2 | dy==3 | dy==4);

betas_over_days = nan(10, 4);
betas_over_days(1:sum(dy==1), 1) = betas_PE(dy==1, 3)';
betas_over_days(1:sum(dy==2), 2) = betas_PE(dy==2, 3)';
betas_over_days(1:sum(dy==3), 3) = betas_PE(dy==3, 3)';
betas_over_days(1:sum(dy==4), 4) = betas_PE(dy==4, 3)';

[a,b,c,d]=ttest2(allBetas1(:, 3), allBetas234(:, 3))
% 2.12, p = 0.0524

[a,b,c,d]=ttest(betas_PE(dy==2, 3)) % 2.83, df = 4, p = 0.047
[a,b,c,d]=ttest(betas_PE(dy==3, 3)) % 1.84, df = 4, p = .14
[a,b,c,d]=ttest(betas_PE(dy==4, 3)) % 3.23, df = 2, p = 0.08

subplot(5, 3, 7)
figInfo = {'Grid-Like Signal', 'Days', 'Betas [a.u.]', [], [-0.1, 0.15], [1 2 3 4], ...
           {'D1', 'D2', 'D3', 'D4'}, [], {}, 24, [], []}; 
[curr] = plotBar(betas_over_days, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, confInt, 0);

subplot(5, 3, 8)
figInfo = {'Day 1', 'Symmetries', 'Betas [a.u.]', [], [-0.1, 0.15], [1 2 3 4 5], ...
           {'4', '5', '6', '7', '8'}, [], {}, 24, [], []}; 
[curr] = plotBar(allBetas1, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, confInt, 0);

subplot(5, 3, 9)
figInfo = {'Day 2-4', 'Symmetries', 'Betas [a.u.]', [], [-0.1, 0.15], [1 2 3 4 5], ...
           {'4', '5', '6', '7', '8'}, [], {}, 24, [], []}; 
[curr] = plotBar(allBetas234, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, confInt, 0);
x = print_ttest(allBetas234(:, 3));
% 0.0010    4.3365   12.0000

subplot(5, 2, 7)
figInfo = {'Day 1', 'Angles', 'z-score [a.u.]', [], [-0.2, 0.2], [1 3 5 7 9 11], {'\omega', '\omega + 60', '\omega + 120', '\omega + 180', '\omega + 240', '\omega + 300'}, [], {}, 24, [], []};
[curr] = plotBar(allSignal1', figInfo, auxStruct.pinkVec, auxStruct.grayVec, curr, confInt, 0);

subplot(5, 2, 8)
figInfo = {'Day 2-4', 'Angles', 'z-score [a.u.]', [], [-0.2, 0.2], [1 3 5 7 9 11], {'\omega', '\omega + 60', '\omega + 120', '\omega + 180', '\omega + 240', '\omega + 300'}, [], {}, 24, [], []};
[curr] = plotBar(allSignal234', figInfo, auxStruct.pinkVec, auxStruct.grayVec, curr, confInt, 0);



%% post-revision

%% grid code x anatomy 

curr = figure(67)

betas = squeeze(nanmean(out_dat1.RegressionResults_glmfit.Test(1, :, 3, :), 4));

% compute based on paxinos / saleem 
ap_vals = auxStruct.chanSessions(:, 7);

px_brain_region = auxStruct.chanSessions(:, 13);
sl_brain_region = auxStruct.chanSessions(:, 14);

px_brain_region(setdiff(1:497, auxStruct.relRegs{4})) = nan;
sl_brain_region(setdiff(1:497, auxStruct.relRegs{4})) = nan;

px_medial  = px_brain_region == 24 | px_brain_region == 32 | px_brain_region == 10 | px_brain_region == 25;

px_ventral = px_brain_region == 14 | px_brain_region == 13;

% both atlases have identical splits here
print_ttest(betas(px_medial))
% 0.0001    4.3283   46.0000

print_ttest(betas(px_ventral))
% 0.0020    3.3581   32.0000

print_ttest(betas(px_medial), 2, betas(px_ventral))
% 0.8188   -0.2298   78.0000

% compute effects for each region
indiv_reg = [24, 32, 10, 25, 14, 13];

px_nums = [];
indx = 1;
for i = indiv_reg
px_nums(indx, :) = [i, sum(px_brain_region == i)];
indx = indx + 1;
end

sl_nums = [];
indx = 1;
for i = indiv_reg
sl_nums(indx, :) = [i, sum(sl_brain_region == i)];
indx = indx + 1;
end

% looking at individual regions for saleem
saleem_effects = []; % for saleem significant 10 / 14 / trending 13
indx = 1;
% looking at individual regions for paxinos 2024
for i = indiv_reg
    o = print_ttest(betas(sl_brain_region == i));
    saleem_effects(indx, :) = [i, o];
    indx = indx + 1;
end
%    24.0000    0.0050  127.3543    1.0000
%    32.0000    0.7772    0.3027    4.0000
%    10.0000    0.0002    4.1623   37.0000
%    25.0000    0.5818    0.7711    1.0000
%    14.0000    0.0100    2.9213   16.0000
%    13.0000    0.0792    1.8831   15.0000

paxinos_effects = [];
indx = 1;
% looking at individual regions for paxinos 2024
for i = indiv_reg
    o = print_ttest(betas(px_brain_region == i));
    paxinos_effects(indx, :) = [i, o];
    indx = indx + 1;
end
%    24.0000    0.0050  127.3543    1.0000
%    32.0000    0.0247    2.5999   11.0000
%    10.0000    0.0000    4.9337   26.0000
%    25.0000    0.4932   -0.7389    5.0000
%    14.0000    0.0002    4.3544   23.0000
%    13.0000    0.7567    0.3206    8.0000

indiv_reg = [32, 10, 14, 13]; % these have the largest number of channels across both atlases

joined_anat_betas = nan(60, 4);
indx = 1;
for i = indiv_reg
    joined_anat_betas(1:length(betas(px_brain_region == i)), indx) = betas(px_brain_region == i);
    indx = indx + 1;
end

subplot(5, 9, 40:41)
figInfo = {{'Paxinos et al. (2024)'}, 'Brodmann Area', 'Betas [a.u.]', [], [-0.25, 0.35], [1 2 3, 4], {'32', '10', '14', '13'}, [], {}, 24, [], []};
[curr] = plotBar(joined_anat_betas, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, [0, 0], 1, 0);

greenVec = [84/255, 161/255, 140/255];
joined_anat_betas = nan(60, 4);
indx = 1;
for i = indiv_reg
    joined_anat_betas(1:length(betas(sl_brain_region == i)), indx) = betas(sl_brain_region == i);
    indx = indx + 1;
end
subplot(5, 9, 42:43)
figInfo = {{'Saleem & Logothetis (2012)'}, 'Brodmann Area', '', [], [-0.25, 0.35], [1 2 3, 4], {'32', '10', '14', '13'}, [], {}, 24, [], []};
[curr] = plotBar(joined_anat_betas, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, [0, 0], 1, 0);

rel_betas = betas(auxStruct.relRegs{4});
ap_vals = auxStruct.chanSessions(auxStruct.relRegs{4}, 7);
site_m = find(ap_vals <= 35.5);
site_a = find(ap_vals > 35.5);

m_betas = rel_betas(site_m);
a_betas = rel_betas(site_a);

% no difference across ap 
% 0.0003    3.8844   46.0000
print_ttest(m_betas)

print_ttest(a_betas)
%  0.0003    4.0119   32.0000

% no difference in comparison
print_ttest(m_betas, 2, a_betas)

joined_anat_betas = nan(60, 2);
joined_anat_betas(1:length(m_betas), 1) = m_betas;
joined_anat_betas(1:length(a_betas), 2) = a_betas;

subplot(5, 9, 44:45)
figInfo = {{''}, 'Anterior-Posterior', '', [0, 3], [0, 0.15], [1, 2], {'<=35', '>35'}, [], {}, 24, [], []};
[curr] = plotBar(joined_anat_betas, figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, [0, 0], 0); axis square

%% cohens'd + concatenated
betas_old = sesBetas{4};
betas_new = out_dat2.vmpfc_pe;

[a,b,c,d]=ttest([betas_old; betas_new(:, 3)])
nanmean([betas_old; betas_new(:, 3)])
nanstd([betas_old; betas_new(:, 3)])

num_sessions = ~isnan([betas_old; betas_new(:, 3)]);
num_sessions = sum(num_sessions);

coh_d = d.tstat/sqrt(num_sessions)
% coh d across both datasets

[a,b,c,d]=ttest(betas_old)
num_sessions = ~isnan([betas_old]);
num_sessions = sum(num_sessions);
coh_d = d.tstat/sqrt(num_sessions)

[a,b,c,d]=ttest(betas_new(:, 3))
num_sessions = ~isnan(betas_new(:, 3));
num_sessions = sum(num_sessions);
coh_d = d.tstat/sqrt(num_sessions)

% how many 
unq_sessions = auxStruct.chanSessions(auxStruct.relRegs{4}, 2);
unq_subs     = auxStruct.chanSessions(auxStruct.relRegs{4}, 6);

% significant sessions per subject
uppVMPFC = squeeze(out_dat1.uppThreshold(:, :, 4));
uppVMPFC = uppVMPFC(1:16, :);
percSig = sum(avgBetasPE > uppVMPFC)./size(uppVMPFC, 1);

rel=find(avgBetasPE(:, 3) > uppVMPFC(:, 3))

unq_s = unique(unq_sessions);
unq_s(rel) % cover boths (1-3 an 0, 4 an 1);
 
end