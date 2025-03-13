function out = figure4_GridNeurons_supplements()

%% load data structs
%cd(go into data folder) 
rel_files = {'unitInfo', 'unitPhaseLocked', 'averagedFRInformation', 'auxStruct', 'gridPhaseBetas', 'AuxiliaryCells', 'averagedPhaseInformation', 'cpdStats', 'stored_betas_roi', 'gridPhaseBetas_sup', 'cpd_inf_perm_vmpfc', 'empir_map_supp'};
for fl = 1:length(rel_files)
    load(rel_files{fl});
end

%% supplements for this figure 

curr = figure(7)
% ----------------------------------
% t-map of phase by phase 
avgedFR_pop = squeeze(nanmean(nanmean(nanmean(unitPhaseLocked, 2), 1), 4))';

curr=figure(666)
subplot(6, 2, [1 3])
tmap = [];
tmap_p = [];
for i = 1:10
    for j = 1:10
        [a,b,c,d]=ttest(avgedFR_pop(:, i), avgedFR_pop(:, j))
        tmap(i, j) = d.tstat;
        tmap_p(i, j) = b;
    end
end

imagesc(triu(tmap, 0)'); 
hcb = colorbar;
hcb.Title
hcb.Title.String = "T-value";
colormap viridis

figElements(curr, 'T-value map', 'Phase', 'Phase',... 
        [0.5 10.5], [0.5 10.5],  [0.5, floor(length(auxStruct.radBounds)/2) - 0.5, length(auxStruct.radBounds)- 0.5], {'0', '\pi', '2\pi'}, [0.5, floor(length(auxStruct.radBounds)/2)-0.5, length(auxStruct.radBounds)-0.5], {'0', '\pi', '2\pi'}, 24, [], []);  

triu(tmap_p)' < (0.05 / 40) % bonferroni correction for 40 perfomed tests on diagonal 

%    0   1   1   1   1   1   1   1   1   1
%    0   0   1   1   1   1   1   1   1   1
%    0   0   0   1   1   1   1   1   1   1
%    1   0   0   0   1   1   1   1   1   1
%    0   0   0   0   0   1   1   1   1   1
%    0   0   0   0   0   0   1   1   1   1
%    0   0   0   0   0   0   0   1   1   1
%    0   0   0   1   1   1   1   0   1   1
%    0   0   1   1   1   1   1   0   0   1
%    0   0   0   1   1   0   0   0   0   0

triu(tmap_p)' * 40

% ----------------------------------
% expand supplements with this figure - preferred / non-prefered phase of cell

% opposite modulation
i = 126;
sample_cells = squeeze(nanmean(nanmean(unitPhaseLocked, 1), 2));
subplot(6, 6, 4:5)
imagesc(imgaussfilt(sample_cells(:, :, i)*1000, 1.0)); hcb = colorbar;
hcb.Title
hcb.Title.String = "Firing Rate [Hz]";
colormap viridis
figElements(curr, 'Cell Example', '', 'Phase',... 
                [1 80], [], [30 50 80], {'-500','Choice ON'}, [1, floor(length(auxStruct.radBounds)/2), length(auxStruct.radBounds)], {'0', '\pi', '2\pi'}, 24, 1, []);  
subplot(6, 6, 6)
sample_cells_trials = squeeze(nanmean(unitPhaseLocked(:, :, 5, :, i), 1));
% remove nans 
sample_cells_trials = removenans(sample_cells_trials, 1);
a=plotmse(sample_cells_trials(:, 1:85) * 1000, [230/255, 159/255, 1/255], [0, 1])
sample_cells_trials = squeeze(nanmean(unitPhaseLocked(:, :, 1, :, i), 1));
sample_cells_trials = removenans(sample_cells_trials, 1);
b=plotmse(sample_cells_trials(:, 1:85) * 1000, [86/255, 180/255, 233/255], [0, 1])
figElements(curr, '', '', 'Firing Rate [Hz]', [1 80], [], [30 50 80], {'-500', 'Choice ON'}, [], {}, 24, 1, []);  
legend([a, b], {'Preferred Phase', 'Non-Preferred Phase'}, 'location', 'best')

% no modulation
i = 89;
subplot(6, 6, 10:11)
imagesc(imgaussfilt(sample_cells(:, :, i)*1000, 1.0)); hcb = colorbar;
hcb.Title
hcb.Title.String = "Firing Rate [Hz]";
colormap viridis
figElements(curr, '', 'Time', 'Phase',... 
                [1 80], [], [30 50 80], {'-500'}, [1, floor(length(auxStruct.radBounds)/2), length(auxStruct.radBounds)], {'0', '\pi', '2\pi'}, 24, [], []);  

subplot(6, 6, 12)
sample_cells_trials = squeeze(nanmean(unitPhaseLocked(:, :, 5, :, i), 1));
% remove nans 
sample_cells_trials = removenans(sample_cells_trials, 1);
a=plotmse(sample_cells_trials(:, 1:85) * 1000, [230/255, 159/255, 1/255], [0, 1])
sample_cells_trials = squeeze(nanmean(unitPhaseLocked(:, :, 10, :, i), 1));
sample_cells_trials = removenans(sample_cells_trials, 1);
b=plotmse(sample_cells_trials(:, 1:85) * 1000, [86/255, 180/255, 233/255], [0, 1])
figElements(curr, '', 'Time [msec]', 'Firing Rate [Hz]', [1 80], [], [30 80], {'-500', 'Choice ON'}, [], {}, 24, [], []);  


% ---
% cpd across individual phases


% ---
% extract the stats for individual phases across cells
permM   = squeeze(nanmean(nanmean(cpdStats), 4));
permSTD = squeeze(nanstd(nanmean(cpdStats), [], 4));

zscores = (squeeze(nanmean(cpdStats(:, :, :, 1))) - permM)./permSTD;
normcdf(-zscores) % pvalue sig at 4th bin

% ---
% plot this for all cells across phase
subplot(6, 3, 8)
cpdSig = squeeze(cpdStats(:, :, :, 1));
cpdSig = squeeze(cpdSig(~isnan(cpdSig(:, 1, 1)), 1, :));
blueVec = auxStruct.blueVec;
b = bar(nanmean(squeeze(cpdSig))*100, 'edgeColor', 'none'); hold on; box off 
b.FaceColor = 'flat';
for j = 1:4
    b.CData(j, :) = blueVec;
end
err = nanstd(squeeze(cpdSig), [], 1)./sqrt(sum(~isnan(squeeze(cpdStats(:, 2, :, 1))), 1));
errorbar(1:4, nanmean(squeeze(cpdSig))*100, err*100,  '.', 'CapSize', 0, 'linewidth', 2, 'color', [155/255, 155/255, 155/255])
figElements(curr, 'VMPFC', 'Seed Window', 'CPD [%]', [], [], [1 2 3 4], {'1', '2', '3', '4'}, [], {}, 24, [], [])

% ---
% temporal alignment; look what is the timeperiod which yields the temporal overlap for the cells
subplot(6, 2, 7)
allModifsAvg = nan(680, 1);
for i = find(auxStruct.brain_region_cells == 5)
    
    allModifsPh = [];
    for j = 1:14
        cx = modifVals{i, j};
        allModifsPh(j, :) = cx;
    end
    allModifsAvg(i) = nanmean(nanmean(allModifsPh));
end

histogram((allModifsAvg + 500) - 800, 'facecolor', auxStruct.grayVec); box off
figElements(curr, '', 'Temporal offset [msec before Choice ON]', 'Frequency',... 
                [], [], [], {}, [], {}, 24, [], []);  

% plot firing rate over time 
% plot individual cells FR
overall_fr = squeeze(nanmean(nanmean(averagedFRInformation(:, :, :, find(auxStruct.brain_region_cells == 5)), 1), 2)) * 1000;

subplot(6, 2, 8)
avgPhase_pl = squeeze(circ_mean(permute(averagedPhaseInformation(:, :, :), [2, 1, 3])));
avgPhase_pl = squeeze(avgPhase_pl(~sum(avgPhase_pl, 2)==0, :));

plotmse((overall_fr(176:476, :) - nanmean(overall_fr, 1))', blueVec, [0, 1])
plot(normalize_bound(circ_mean(avgPhase_pl(:, 176:476))', -0.2, 0.2), 'k', 'linewidth', 1.5); hold on
figElements(curr, 'Population Firing rate', 'Time [msec]', 'Firing Rate Change [Hz]',... 
                [1 270], [], [], {}, [], {}, 24, [], []);  

% ---
% look whether griddy cells are related to cells encoding chosen val
relCPD    = squeeze(cpdStats(find(auxStruct.brain_region_cells == 5), :, :, 1));
relGRID   = squeeze(stored_betas_roi(find(auxStruct.brain_region_cells == 5), :, 2)); % 
  
corr(squeeze(relCPD(:, 1, :)), relGRID(:, 3))
% do also t-test according to cpd based on gridness
high_grid = find(relGRID(:, 3) > nanmedian(relGRID(:, 3)));
low_grid  = find(relGRID(:, 3) <= nanmedian(relGRID(:, 3))); 
relCPD_info = squeeze(relCPD(:, 1, :)); % also not true in t-test
[a,b,c,d]=ttest2(relCPD_info(high_grid), relCPD_info(low_grid))

subplot(6, 3, 9)
scatterfit(squeeze(relCPD(:, 1, 4)), relGRID(:, 3), 1)
figElements(curr, '', 'Chosen Value CPD [%]', 'Hexadirectional Modulation Beta [a.u.]', [], [], [], {}, [], {}, 24, [], []);

% -------
% we compare the other clusters
cell_cpd = cpd_inf_perm(:, :, 1);
relCPD  = squeeze(cell_cpd);
relGRID = permute(squeeze(nanmean(gridPhaseBetas(:, 1, :, :, find(auxStruct.brain_region_cells == 5)), 1)), [3, 1, 2]); relGRID = relGRID(:, :, 3);

only_full = find(~any(isnan(relCPD), 2));

subplot(6, 3, 13)
scatterfit(relCPD(only_full, 245), relGRID(only_full, 228), 1)
figElements(curr, '1st Cluster', 'Chosen Value CPD [%]', 'Hexadirectional Modulation Beta [a.u.]', [], [], [], {}, [], {}, 24, [], []);
subplot(6, 3, 14)
scatterfit(relCPD(only_full, 343), relGRID(only_full, 228), 1)
figElements(curr, '2nd Cluster', 'Chosen Value CPD [%]', 'Hexadirectional Modulation Beta [a.u.]', [], [], [], {}, [], {}, 24, [], []);
subplot(6, 3, 15)
scatterfit(relCPD(only_full, 392), relGRID(only_full, 228), 1)
figElements(curr, '3rd Cluster', 'Chosen Value CPD [%]', 'Hexadirectional Modulation Beta [a.u.]', [], [], [], {}, [], {}, 24, [], []);

% ---
% hexadir analysis in a phase dependent way using the bins 

subplot(6, 3, 7)
tmp = permute(squeeze(nanmean(gridPhaseBetas_sup(:, :, :, 1:sum(auxStruct.brain_region_cells == 5)), 1)), [3, 2, 1]);
[a,b,c,d]= ttest(tmp);
tstats = squeeze(d.tstat)

%     0.5255    0.6131    0.0415    0.3229    0.3702    0.8079    0.1045    0.2295    0.6310    0.3601
%     0.9034    0.1871    0.4156    0.2356    0.4430    0.8784    0.6117    0.5820    0.1762    0.7194
%     0.6363    0.4385    0.2880    0.8289    0.3541    0.0249    0.2384    0.8511    0.5101    0.7667
%     0.3773    0.8118    0.8628    0.3866    0.0511    0.1418    0.1790    0.6455    0.9220    0.1776
%     0.7721    0.7535    0.1596    0.1568    0.9824    0.5701    0.9518    0.1363    0.9253    0.8543

%     0.6362    0.5067    2.0556    0.9917    0.8986    0.2436   -1.6329   -1.2063   -0.4812   -0.9178
%    -0.1215   -1.3249   -0.8163   -1.1907    0.7691    0.1532    0.5086   -0.5515   -1.3585    0.3599
%    -0.4738    0.7766    1.0661    0.2165    0.9295    2.2640    1.1834    0.1880   -0.6602    0.2972
%     0.8853   -0.2385    0.1730    0.8682    1.9651    1.4763    1.3498    0.4609    0.0981   -1.3541
%     0.2902   -0.3146   -1.4130   -1.4227   -0.0221   -0.5691   -0.0605   -1.4973   -0.0939   -0.1839

tstatsStatic = tstats;

% add extra val for visual purposes to keep figure wrapped

filler = tstatsStatic(:, 1);
tstatsStatic(:, 11) = filler;
plot(squeeze(tstatsStatic([1 2 4 5],  :))', 'color', auxStruct.grayVec, 'linewidth', 1); hold on; box off
plot(squeeze(tstatsStatic([3], :))', 'color', auxStruct.greenVec, 'linewidth', 3.5); hold on %

plot(1:length(auxStruct.radBounds), repmat(-1.975, [length(auxStruct.radBounds), 1]), '--k');
plot(1:length(auxStruct.radBounds), repmat(1.975, [length(auxStruct.radBounds), 1]), '--k');

legend('4-fold', '5-fold', '7-fold', '8-fold', '6-fold'); legend off
figElements(curr, '', 'Phase', 'Hexadirectional Modulation [T-value]',...
    [1 11], [-3 3], [1 6, 11], {'0', '\pi', '2\pi'}, [], {}, 24, [], []);


%% post-revision analyses


%% evaluate chosen value coding results 

% extract betas
out_betas = squeeze(nanmean(gridPhaseBetas(:, 1, 1:476, :, find(auxStruct.brain_region_cells == 5))));    
bet_v = out_betas(:, :, :);
mean_beta = squeeze(nanmean(bet_v(auxStruct.grid_peak, :, :), 1))';
[~,rnkd_vals] = sort(mean_beta(:, 3) - sum(abs(mean_beta(:, [1, 2, 4, 5])), 2), 'descend');

warning off
% split in half
top_third    = rnkd_vals(1:round(auxStruct.num_vmpfc_cells/3));
bottom_third = rnkd_vals(auxStruct.num_vmpfc_cells - length(top_third) + 1:end);

%% fourier results

% preallocate
unshifted_fourier = nan(5, 5, sum(auxStruct.brain_region_cells == 5), auxStruct.nPerm + 1);

for np = 1:(auxStruct.nPerm + 1)
    for i = 1:sum(auxStruct.brain_region_cells == 5)
        curr_map = data_file(:, :, i, np);
        curr_map(isnan(curr_map)) = 0; % zero out nans otherwise cannot be computed
        curr_map = reshape(curr_map, [5, 5]);  
        fft1 = abs(fft2(curr_map)); % take abs magnitude 
        unshifted_fourier(:, :, i, np) = fft1;
    end
end

foo = reshape(unshifted_fourier, [25, 1, sum(auxStruct.brain_region_cells == 5), auxStruct.nPerm + 1]);
% compute freqs along x and y coord 
% freqs are split based on median because this captures relationship between non-neighbouring states (high) and neighbouring states (low)
% this is specific for this space size and doesn't work for larger spaces
n_states = 5;
fx = (0:n_states-1) / n_states;  
fx(fx > 0.5) = fx(fx > 0.5) - 1;
radius_unshifted = sqrt(fx.^2 + (fx.^2)'); 
radius = radius_unshifted;

high_freq_mask = radius >= median(radius(:)); % above .4 
low_freq_mask  = radius < median(radius(:)); 
low_freq_mask(1, 1) = 0; % dc offset not included in comparison

far   = nanmean(foo(find(high_freq_mask(:)), :, :, 1), 1);
close = nanmean(foo(find(low_freq_mask(:)), :, :, 1), 1);
dc_offset = squeeze(foo(1, :, :, 1));

ratio = squeeze((close - far)); 
print_ttest(ratio(top_third))
% 0.0222    2.3573   52.0000
print_ttest(ratio(bottom_third))
% 0.8403    0.2026   52.0000

print_ttest(ratio(top_third), 2, ratio(bottom_third))
% 0.1323    1.5169  104.0000

rel_fourier = squeeze(ratio)'; 
rel_betas = mean_beta(:, 3);  % hex signal
tmp_c = find(~isnan(rel_betas));
rel_betas   = rel_betas(tmp_c);
rel_fourier = rel_fourier(tmp_c);

tmp_b = find(~isnan(rel_fourier));
rel_betas   = rel_betas(tmp_b);
rel_fourier = rel_fourier(tmp_b);

curr = figure(6)

% lots of outliers so have to use spearman to minimize impact
[a,b] = corr(rel_betas, rel_fourier', 'type', 'spearman')

subplot(2, 3, 2)
scatterfit(rel_betas, rel_fourier', 1, auxStruct.greenVec, 'spearman')
figElements(curr, '', 'Hexadirectional Modulation Beta [a.u.]', 'Fourier Magnitude Difference [a.u.]', [], [], [], {}, [], {}, 24, [], []);

orangeCol = [0.8, 0.6, 0.0];

curr = figure(8)
subplot(2, 3, 3)
mean_vals_all = nan(160, 2);
ratio = squeeze(close - far);
mean_vals_all(1:size(top_third, 1), 1)    = ratio(top_third); 
mean_vals_all(1:size(bottom_third, 1), 2) = ratio(bottom_third);

figInfo = {'', 'Grid-like Code', 'Fourier Frequency Difference [a.u.]', ... 
                 [], [], [1, 2], {'High', 'Low'}, [], {}, 24, [], []};
plotBar(mean_vals_all, figInfo, auxStruct.greenVec, orangeCol, curr, [0, 1], 0); 


% compare all hypotheses
num_hypotheses = 17;

% define tmp ev 
tmp_ev = flipud((1:5) .* (5:-1:1)');
tmp    = tmp_ev;

all_hypotheses = nan(5, 5, num_hypotheses);
value_differ = nan(auxStruct.num_vmpfc_cells, 90, num_hypotheses);

% partial checkerboard
for k  = 1:num_hypotheses
    c_ind = 1;
    for kk = 1:auxStruct.num_vmpfc_cells
        curr_mean_fr = data_file(:, :, kk, 1);
        curr_mean_fr(isnan(curr_mean_fr)) = 0;
        curr_mean_fr = nanmean(curr_mean_fr(:));
    
        % define all four
    
        if k == 1
            % checker
            tmp_check = zeros(5, 5);
            tmp_check(1:2:end) = 1;
            current_avg = mean(tmp_check(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_check = tmp_check .* scaling_factor;
        
            tmp = tmp_check;

        elseif k == 2
            % partial checker
            tmp_partial_check = zeros(5, 5);
            tmp_partial_check(1:2:end) = 1;
            tmp_partial_check(1:2, :) = 0;
            current_avg = mean(tmp_partial_check(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_partial_check = tmp_partial_check .* scaling_factor;
        
            tmp = tmp_partial_check;

        elseif k == 3

            % quadrant checker
            tmp_quadrant_check = zeros(5, 5);
            tmp_quadrant_check(1:2:end) = 1;
            tmp_quadrant_check(tmp_ev >= nanmedian(tmp_ev(:))==0) = 0;
            current_avg = mean(tmp_quadrant_check(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_quadrant_check = tmp_quadrant_check .* scaling_factor;
        
            tmp = tmp_quadrant_check;    

        elseif k == 4
            % hex
            tmp_hex = zeros(5, 5);
            tmp_hex(1, [2, 4]) = 1;
            tmp_hex(3, [1, 3, 5]) = 1;
            tmp_hex(5, [2, 4]) = 1;
            current_avg = mean(tmp_hex(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_hex = tmp_hex .* scaling_factor;
        
            tmp = tmp_hex;          

        elseif k == 5
            % partial hex 
            tmp_partial_hex = zeros(5, 5);
            tmp_partial_hex(3, [1, 3, 5]) = 1;
            tmp_partial_hex(5, [2, 4]) = 1;
            current_avg = mean(tmp_partial_hex(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_partial_hex = tmp_partial_hex .* scaling_factor;
        
            tmp = tmp_partial_hex; 

        elseif k == 6
            % quad hex 
            tmp_quad_hex = zeros(5, 5);
            tmp_quad_hex(1, [2, 4]) = 1;
            tmp_quad_hex(3, [1, 3, 5]) = 1;
            tmp_quad_hex(5, [2, 4]) = 1;
            tmp_quad_hex(tmp_ev >= nanmedian(tmp_ev(:))==0) = 0;
            current_avg = mean(tmp_quad_hex(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_quad_hex = tmp_quad_hex .* scaling_factor;
        
            tmp = tmp_quad_hex; 

        elseif k == 7

            % sparse hex
            tmp_sparse_hex = zeros(5, 5);
            tmp_sparse_hex(1:2:end) = 1;
            tmp_sparse_hex([2, 4, 5], :) = 0;
            tmp_sparse_hex(1, [1, 5]) = 0;
            tmp_sparse_hex(3, 3) = 0;
            current_avg = mean(tmp_sparse_hex(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_sparse_hex = tmp_sparse_hex .* scaling_factor;
        
            tmp = tmp_sparse_hex;            


        elseif k == 8
            % mag 
            tmp_mg = repmat((1:5)', [1, 5]);
            current_avg = mean(tmp_mg(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_mg = tmp_mg .* scaling_factor;

            tmp = tmp_mg';

        elseif k == 9
            % partial_mag
            tmp_mg = repmat((1:5)', [1, 5]);
            tmp_mg(1:2, :) = 0;
            current_avg = mean(tmp_mg(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_mg = tmp_mg .* scaling_factor;

            tmp = tmp_mg';
  
        elseif k == 10
            % quad mag 
            tmp_mg = repmat((1:5)', [1, 5]);
            tmp_mg(tmp_ev >= nanmedian(tmp_ev(:))==0) = 0;
            current_avg = mean(tmp_mg(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_mg = tmp_mg .* scaling_factor;

            tmp = tmp_mg';

        elseif k == 11
            % tmp ev 
            tmp_ev = flipud((1:5) .* (5:-1:1)');
            current_avg = mean(tmp_ev(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_ev = tmp_ev .* scaling_factor;
            
            tmp = tmp_ev;
        elseif k == 12
            
            % tmp partial ev 
            tmp_ev = flipud((1:5) .* (5:-1:1)');
            tmp_ev(1:2, :) = 0;
            current_avg = mean(tmp_ev(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_ev = tmp_ev .* scaling_factor;
            
            tmp = tmp_ev;
        elseif k == 13
            % tmp quad ev 
            tmp_ev = flipud((1:5) .* (5:-1:1)');
            tmp_ev(tmp_ev >= nanmedian(tmp_ev(:))==0) = 0;
            current_avg = mean(tmp_ev(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_ev = tmp_ev .* scaling_factor;
            
            tmp = tmp_ev;

        elseif k == 14

            % high value
            tmp_high_val= zeros(5, 5);
            tmp_high_val(4, 4:5) = 1;
            tmp_high_val(5, 4:5) = 1;
            current_avg = mean(tmp_high_val(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_high_val = tmp_high_val .* scaling_factor;
        
            tmp = tmp_high_val;

        elseif k == 15

            % quadratic code
            tmp_quadrat = zeros(5, 5);
            tmp_quadrat(1, :) = 2;
            tmp_quadrat(2, :) = 1;
            tmp_quadrat(3, :) = 0;
            tmp_quadrat(4, :) = 1;
            tmp_quadrat(5, :) = 2;
            current_avg = mean(tmp_quadrat(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_quadrat = tmp_quadrat .* scaling_factor;
        
            tmp = tmp_quadrat;

        elseif k == 16
            % quadratic partial code
            tmp_quadrat_p = zeros(5, 5);
            tmp_quadrat_p(1, :) = 0;
            tmp_quadrat_p(2, :) = 0;
            tmp_quadrat_p(3, :) = 0;
            tmp_quadrat_p(4, :) = 1;
            tmp_quadrat_p(5, :) = 2;
            current_avg = mean(tmp_quadrat_p(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_quadrat_p = tmp_quadrat_p .* scaling_factor;
        
            tmp = tmp_quadrat_p;

        elseif k == 17
            % quadratic quad code
            tmp_quadrat_q = zeros(5, 5);
            tmp_quadrat_q(1, :) = 2;
            tmp_quadrat_q(2, :) = 1;
            tmp_quadrat_q(3, :) = 0;
            tmp_quadrat_q(4, :) = 1;
            tmp_quadrat_q(5, :) = 2;
            tmp_quadrat_q(tmp_ev >= nanmedian(tmp_ev(:))==0) = 0;
            current_avg = mean(tmp_quadrat_q(:));
            scaling_factor = curr_mean_fr ./ current_avg;
            tmp_quadrat_q = tmp_quadrat_q .* scaling_factor;
        
            tmp = tmp_quadrat_q;            
                  
        end

        all_hypotheses(:, :, k) = tmp;

        for i = 1:1:180
        fft1 = fft2(imrotate(tmp, i, 'bilinear', 'crop'));
        mag1 = abs(fft1);

        high_freq_mask = radius >= median(radius(:)); 
        low_freq_mask  = radius < median(radius(:)); 
        low_freq_mask(1, 1) = 0;

        far   = nanmean(mag1(find(high_freq_mask(:))), 1);
        close = nanmean(mag1(find(low_freq_mask(:))), 1);

        differ_theor = (close-far); 
        value_differ(c_ind, i, k) = differ_theor; 
        end
        c_ind = c_ind + 1;
    end

    k
end

vmpfc_cells = find(auxStruct.brain_region_cells == 5);

target_cells = 1:length(vmpfc_cells);
foo = reshape(squeeze(unshifted_fourier(:, :, :, 1)), [25, length(vmpfc_cells)]);
far   = nanmean(foo(find(high_freq_mask(:)), target_cells), 1);
close = nanmean(foo(find(low_freq_mask(:)), target_cells), 1);
target_fourier_vals = (close-far);

[a, mx] = min((squeeze(nanmean(value_differ(target_cells, :, :), 1))-nanmean(target_fourier_vals)).^2);
[~, my] = min(a);

% best fitting template
tmp = zeros(5, 5);
tmp(1:2:end) = 1;

subplot(2, 2, 1)
imagesc(imrotate(tmp, mx(my), 'bilinear', 'crop')); colormap viridis
figElements(curr, 'Best Fitting Template', 'Magnitude', 'Probability', [], [], [1, 5], {'Min', 'Max'}, [1, 5], {'Min', 'Max'}, 8, [], []);
% 
subplot(2, 2, 2)
fft1 = fft2(imrotate(tmp,  mx(my), 'bilinear', 'crop'));
mag1 = abs(fft1);
mag1(1, 1) = nan;
imagesc(mag1)
title('theoretical coefficients')

%differences in fr 
top_fr = squeeze(nanmean(nanmean(nanmean(averagedFRInformation(:, :, auxStruct.grid_peak, vmpfc_cells(top_third)), 3), 2), 1));
bt_fr  = squeeze(nanmean(nanmean(nanmean(averagedFRInformation(:, :, auxStruct.grid_peak, vmpfc_cells(bottom_third)), 3), 2), 1));

% not significantly different
[a,b,c,d] = ttest2(top_fr, bt_fr)

% pca results 
rel_map = data_file(:, :, :, 1);
rel_map = reshape(rel_map, [25, auxStruct.num_vmpfc_cells]);
rel_map(isnan(rel_map)) = 0;

% remove global mean effect
rel_map = rel_map - nanmean(rel_map, 2);

differ_all = [];
for i = 1:(auxStruct.nPerm + 1)
    if i == 1
        cell_a = top_third;
        cell_b =  bottom_third;
    elseif i > 1
        cell_all = shuffle([top_third;bottom_third]);
        cell_a = cell_all(1:53); % take equal amt
        cell_b = cell_all(54:end);
    end

    rel_vals = rel_map(:, cell_a);
    [a,b,c, ~, e]=pca(rel_vals');

    rel_vals = rel_map(:, cell_b);
    [a,b,c_a,~, e_a]=pca(rel_vals');
    
    differ = cumsum(e_a(1:25)) - cumsum(e(1:25));
    
    differ_all(:, i) = differ;
i
end

curr = figure(6)
orangeCol = [0.8, 0.6, 0.0];
subplot(2, 3, 6)
cell_a  = top_third;
cell_b =  bottom_third;
rel_vals = rel_map(:, cell_a); 
rel_vals(isnan(rel_vals)) = 0;
[a,b,c, ~, e]=pca(rel_vals');
rel_vals = rel_map(:, cell_b); 
rel_vals(isnan(rel_vals)) = 0;
[a,b,c_a,~, e_a]=pca(rel_vals');
a=plot(cumsum(e), 'color', auxStruct.greenVec, 'linewidth', 2); hold on; box off
b=plot(cumsum(e_a), 'color', orangeCol, 'linewidth', 2)
sig_vec = nan(1, 25);
sig_vec(find(differ_all(:, 1)>prctile(differ_all(:, 2:end), 95, 2) == 1)) = 75;
c = plot(sig_vec, 'color', 'k', 'linewidth', 2)
legend([a, b], {'High', 'Low'})
figElements(curr, '', 'Number of Principal Components', '% Variance Explained', [], [], [], {}, [], {}, 24, [], []);


[a]=(sum(differ_all(1:25, 1)) - nanmean(sum(differ_all(1:25, 2:end))))./nanstd(sum(differ_all(1:25, 2:end)));
p = (1 - normcdf(abs(a)))


%% cells exhibit more periodic activity as a function of distance between options  

cell_indx     = 1;
avg_vector_fr = nan(7, 100, length(vmpfc_cells), 2);
rt_vals       = nan(7, length(vmpfc_cells));

for u = vmpfc_cells
    rel_fr  = squeeze(averagedFRInformation(:, :, :, u));

    left_prob  = discretize(AuxiliaryCells(u).Left_prob, 5);
    right_prob = discretize(AuxiliaryCells(u).Right_prob, 5);
    left_mag   = discretize(AuxiliaryCells(u).Left_pay, 5);
    right_mag  = discretize(AuxiliaryCells(u).Right_pay, 5);

    trial_len = round(sqrt(sum([left_mag - right_mag, left_prob-right_prob].^2, 2))); % euclidean
    
     foo = state_mean(trial_len, squeeze(nanmean(rel_fr(:, 1:length(trial_len), 186:285), 1))); % 100 msec window centred on grid code ~ 186:285 
     avg_vector_fr(foo(:, 1, 1) + 1, :, cell_indx) = squeeze(foo(:, 2, :));
     tr_nums(foo(:, 1, 1) + 1, cell_indx) = squeeze(foo(:, 3, 1));    

    % reaction time 
    trialRT = [];
    for j = 1:length(AuxiliaryCells(u).Codes)
        foo = AuxiliaryCells(u).Codes{j};
        trialRT(j) = foo(foo(:, 1) == 4, 2) - foo(foo(:, 1) == 62, 2);
    end    
    foo = state_mean(trial_len, trialRT');
    rt_vals(foo(:, 1, 1) + 1, cell_indx) = squeeze(foo(:, 2, 1));
    cell_indx = cell_indx + 1;

    u
end

% get fourier per cell per distance
estimated_fft      = nan(6, auxStruct.num_vmpfc_cells, 1);
estimated_fft_full = nan(6, auxStruct.num_vmpfc_cells, 50, 2);
half_point = round(size(avg_vector_fr, 2))/2;

for tp = 1
for i = 1:6
    for j = 1:length(vmpfc_cells)
        curr_val = zscore(avg_vector_fr(i, :, j, tp));
        [ou] = fft(curr_val');

        fft1_shifted = fftshift(ou);

        mag_fft1 = abs(fft1_shifted);
        fft_estimate = mag_fft1((half_point + 1):end); 

        [~,~,st] = glmfit(1:half_point, fft_estimate);

        estimated_fft(i, j, tp) = st.beta(2);
        estimated_fft_full(i, j, :, tp) = fft_estimate;
    end
end
end

% compare frequency decay in top vs. bottom 
[a,b,c,d]=ttest(nanmean(estimated_fft(4:6, :, 1)), nanmean(estimated_fft(1:3, :, 1)))
[a,b,c,d]=ttest(nanmean(estimated_fft(4:6, top_third, 1)), nanmean(estimated_fft(1:3, top_third, 1)))
[a,b,c,d]=ttest(nanmean(estimated_fft(4:6, bottom_third, 1)), nanmean(estimated_fft(1:3, bottom_third, 1)))

long_estimates  = nanmean(estimated_fft(4:6, :))';
short_estimates = nanmean(estimated_fft(1:3, :))';

curr = figure(7)
subplot(2, 4, 2)
% sort for example pattern
[~,zo]=sort(nanmean(estimated_fft(4:6, top_third))-nanmean(estimated_fft(1:3, top_third)), 'descend');

plot(squeeze(nanmean(estimated_fft_full(1:3, top_third(zo(4)), :))), 'linewidth', 1, 'color', auxStruct.shortCol); hold on
plot(squeeze(nanmean(estimated_fft_full(4:6, top_third(zo(4)), :))), 'linewidth', 1, 'color', auxStruct.longCol);
scatter(1:25, squeeze(nanmean(estimated_fft_full(1:3, top_third(zo(4)), 1:25))), 'filled', 'markerfacecolor', auxStruct.shortCol); hold on;

scatter(1:25, squeeze(nanmean(estimated_fft_full(4:6, top_third(zo(4)), 1:25))), 'filled', 'markerfacecolor', auxStruct.longCol); hold on; box off
figElements(curr, '', 'Frequency', 'Normalized Magnitude', [1, 25], [-5 65], [1:10:50], {f(1:10:50)}, [], {}, 24, [], []);

subplot(2, 4, 3)
figInfo = {'Periodicity', '', 'Beta Diff_{Long - Short} [a.u.]', ... 
                 [0, 2], [], [1], {'All Neurons'}, [], {}, 24, [], []};
plotBar([long_estimates-short_estimates], figInfo, auxStruct.greenVec, auxStruct.greenVec, curr, [0, 1], 0); ylim([-0.004, .035]); 

subplot(2, 4, 4)
mean_vals_all = nan(60, 2);
mean_vals = [nanmean(estimated_fft(4:6, bottom_third))'-nanmean(estimated_fft(1:3, bottom_third))'];
mean_vals_all(1:size(mean_vals, 1), 2) = mean_vals;

mean_vals = [nanmean(estimated_fft(4:6, top_third))'-nanmean(estimated_fft(1:3, top_third))'];
mean_vals_all(1:size(mean_vals, 1), 1) = mean_vals;

figInfo = {'Periodicity', 'Grid-like Code', '', ... 
                [-0.5, 3.5], [], [1, 2], {'High', 'Low'}, [], {}, 24, [], []};
plotBar(mean_vals_all, figInfo, auxStruct.greenVec, orangeCol, curr, [0, 1], 0); ylim([-0.004, .035]); 

% difference between groups
print_ttest(nanmean(estimated_fft(4:6, top_third))-nanmean(estimated_fft(1:3, top_third)), 2, nanmean(estimated_fft(4:6, bottom_third))-nanmean(estimated_fft(1:3, bottom_third)))
% 0.0262    2.2560  104.0000

% example indiv cells 
avg_vector_info = {};
indx = 1;
for kk = 1:length(vmpfc_cells)
    
    rel_fr  = squeeze(averagedFRInformation(:, :, :, vmpfc_cells(kk)));
    
    left_prob  = discretize(AuxiliaryCells(vmpfc_cells(kk)).Left_prob, 5);
    right_prob = discretize(AuxiliaryCells(vmpfc_cells(kk)).Right_prob, 5);
    left_mag   = discretize(AuxiliaryCells(vmpfc_cells(kk)).Left_pay, 5);
    right_mag  = discretize(AuxiliaryCells(vmpfc_cells(kk)).Right_pay, 5);
    
    trial_len = round(sqrt(sum([left_mag - right_mag, left_prob-right_prob].^2, 2))); % euclidean
    trial_type = ones(length(left_prob), 2);


    for tr_tp = 1:2
        for tx = unique(trial_len)'
            avg_vector_info{tx+1, indx, tr_tp} = squeeze(nanmean(rel_fr(:, find(trial_len == tx & trial_type(:, tr_tp) == 1), 136:335), 1)); % take longer window to better visualise
        end
    end
    indx = indx + 1;
    kk
end

indx = 1;
for i = top_third([10, 42, 6])' 

    avged_vals = [];
    for k = 1:6
        avged_vals(k, :) = nanmean(avg_vector_info{k, i, 1}, 1); 
    end

    subplot(2, 3, indx + 3)
    long_fr  = nanmean(squeeze(avged_vals(4:6, :)), 1)*1000;
    short_fr = nanmean(squeeze(avged_vals(1:3, :)), 1)*1000;

    rectangle('Position', [85, min([long_fr, short_fr]), 31, max([long_fr, short_fr])], ...
    'FaceColor', [auxStruct.greenVec, 0.1], ...
    'EdgeColor', 'none'); hold on

   b=plot(short_fr, 'linewidth', 3, 'color', auxStruct.shortCol); hold on; box off
   a=plot(long_fr, 'linewidth', 3, 'color', auxStruct.longCol);    
    if i == 34
        legend([a,b], {'Map Dist_{Long}', 'Map Dist_{Short}'}, 'location', 'best')
    end
    figElements(curr, '', 'Time [msec]', 'Firing Rate [Hz]', [], [min([long_fr, short_fr]), max([long_fr, short_fr])], [], {}, [], {}, 24, [], []);
    indx = indx + 1;

end

%% gen fig for lower part

curr = figure(7)
% hypotheses examples
for i = 1:17
names = {'Square Grid', 'High Attr. Value x Square Grid', 'High EV x Square Grid', ...
        'Hexagonal Grid', 'High Attr. Value x Hexagonal Grid', 'High EV x Hexagonal Grid', 'Sparse Hexagonal Grid', ...
        'Linear Attribute Code', 'High Attr. Value x Linear Attribute', 'High EV x Linear Attribute', ...
        'EV Code', 'High Attr. Value x Expected Value', 'High EV x Expected Value', ...
        'Highest Value Code', ...
        'Quadratic Value Code', 'High Attr. Value x Quadratic Value Code', 'High EV x Quadratic Value Code'};
subplot(3, 6, i)
imagesc(all_hypotheses(:, :, i)); colormap viridis; 
figElements(curr, names{i}, 'Magnitude', 'Probability', [], [], [1, 5], {'Min', 'Max'}, [1, 5], {'Min', 'Max'}, 14, [], []);
end


end
