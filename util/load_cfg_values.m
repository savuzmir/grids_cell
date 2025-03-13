function [auxStruct, chanInfo] = load_cfg_values(auxLFP, AuxiliaryCells, baseP)

% set all variables here so they are accessed from one spot 
auxStruct     = struct;
files_to_load = {'all_units_info', 'dayInformation_stimSet', 'channelTrials'};

load('C:\Users\sebas\My Drive\Work\KennerleyLab\Code\Infogathering\cell_paper_revision\master\data\100MsecMT_ResponseMatrixLFP_highThetaB.mat');

for fl = 1:length(files_to_load)
    full_path = fullfile(baseP, [files_to_load{fl}, '.mat']);
    tmp = load(full_path);
    if fl == 1
        auxStruct.(files_to_load{fl}) = tmp.all_units;
    elseif fl == 2
        auxStruct.(files_to_load{fl}) = tmp.both;
    elseif fl == 3
        auxStruct.(files_to_load{fl}) = tmp.channelTrials;
    end
end

auxStruct.nPerm       = 1000;          % number of permutations
auxStruct.crosvals    = 3;             % number of crosvalidation runs
auxStruct.symmetries  = 4:8;           % number of symmetries for hex
auxStruct.num_symm    = 5;
auxStruct.num_sessions_dat1 = 54;

auxStruct.windowSize    = 100;       % size of sliding window
auxStruct.slideWidth    = 10;        % incremented by
auxStruct.epochDur      = 2500;      % duration of full epoch for neural fr raster
auxStruct.preEvent      = 799;       % msec before onset starts
auxStruct.num_channels  = 497; 
auxStruct.num_cells     = 724;
auxStruct.areaCodes     = [1, 2, 3, 5];

% define colors
auxStruct.grayVec            = [148/255, 149/255, 153/255];
auxStruct.greenVec           = [84/255, 161/255, 140/255];
auxStruct.blueVec            = [1/255,  114/255, 178/255];
auxStruct.yellowVec          = [230/255, 159/255, 1/255];
auxStruct.pinkVec            = [243/255, 127/255, 128/255];
auxStruct.areaCol            = [237 228 22;62 181 156;53 45 132; 1 1 1]/255;
longCol  = ([247, 143, 129] - 50)./255;
shortCol = ([127, 128, 200] - 50)./255;
auxStruct.longCol  = longCol;
auxStruct.shortCol = shortCol;
curr_palette = load([baseP, 'RdGy.mat']);
curr_palette = flip(curr_palette.cmap_interp);
auxStruct.curr_palette = curr_palette;

% extract relevant info 
animal  = [];
nTR     = [];
session = [];

for i = 1:auxStruct.num_channels
    nTR(i)                   = size(auxLFP(i).ch.ChEVL, 1);
    animal(i)                = auxLFP(i).ch.sb == 'F';
    session(i)               = chanInfo(i).Aux.session(i);
end

brain_region_channels = [];

for i = 1:length(auxLFP)
    brain_region_channels(i) = auxLFP(i).ch.brain_region;
end

brain_region_cells = [];
for i = 1:length(AuxiliaryCells)
    brain_region_cells(i) = AuxiliaryCells(i).brain_region;
end

% min amt trials for cros-val lfp  analysis in dat1
minTrials = 200;

auxStruct.trialLim               = minTrials;
auxStruct.brain_region_channels  = brain_region_channels;
auxStruct.brain_region_cells     = brain_region_cells;
auxStruct.nTR                    = nTR;
auxStruct.animal                 = animal;
auxStruct.session                = session;

% clear bounds for data 
auxStruct.slSizeStart   =   -800;  % onset 
auxStruct.slSizeEnd     =   800;   % offset 

startBin = round((auxStruct.preEvent + auxStruct.slSizeStart)/auxStruct.slideWidth + 1);
endBin = round((auxStruct.preEvent + auxStruct.slSizeEnd)/auxStruct.slideWidth + 1);

TimePoints = [startBin:endBin];

% bounds for analysis
auxStruct.slSizeStart   =   -300;                                        
auxStruct.slSizeEnd     =   0;                                       

startBin = round((auxStruct.preEvent + auxStruct.slSizeStart)/auxStruct.slideWidth + 1);
endBin   = round((auxStruct.preEvent + auxStruct.slSizeEnd)/auxStruct.slideWidth + 1);

auxStruct.startBin = startBin;
auxStruct.endBin   = endBin;

auxStruct.timepoints_sr      = TimePoints;     % timepoints for regression 
auxStruct.timepoints_sr_perm = 1:121;          % timepoints for permutations
auxStruct.num_tp_dat2        = 300;
auxStruct.permThreshold      = 95;
auxStruct.peakVMPFC          = 63:73;

% define rel values for regions
rgIdx = 1;
for rg = auxStruct.areaCodes
    rel_reg = find(brain_region_channels == rg);
    auxStruct.relRegs{rgIdx}    = intersect(rel_reg, find(nTR > minTrials));
    rgIdx  = rgIdx + 1;
end

% define colors
auxStruct.grayVec            = [148/255, 149/255, 153/255];
auxStruct.greenVec           = [84/255, 161/255, 140/255];
auxStruct.blueVec            = [1/255,  114/255, 178/255];
auxStruct.yellowVec          = [230/255, 159/255, 1/255];
auxStruct.pinkVec            = [243/255, 127/255, 128/255];
% auxStruct.sigVMPFC           = [175   176   232   234   250   251   252   253   399   406   434   435   445   449];
auxStruct.areaCol            = [237 228 22;62 181 156;53 45 132; 1 1 1]/255;


% baseP_Alt = 'C:\Users\sebas\OneDrive - University College London\Sebastijan\KennerleyLab\Data\FIL_update\simultaneous_trials';    
% fixMatFR  = load([baseP_Alt, '\0_raw_FR\FixMatrix.mat']);                    fixMatFR  = fixMatFR.ContainerFR;
% respMatFR = load([baseP_Alt, '\0_raw_FR\ResponseMatrix.mat']);               respMatFR = respMatFR.ContainerFR;
% respMat   = load([baseP_Alt, '\0_raw_LFP_reduced\ResponseMatrixLFP.mat']);   respMat   = respMat.currLFP;
% fixMat    = load([baseP_Alt, '\0_raw_LFP_reduced\FixMatrixLFP.mat']);        fixMat    = fixMat.currLFP;
% rewMat    = load([baseP_Alt, '\0_raw_LFP_reduced\RewardMatrixLFP.mat']);     rewMat    = rewMat.currLFP;

% compile hex regressors for lfp 
sessRT = {};
rewarded = {};
for i = 1:auxStruct.num_channels
    
    % compute rt
    for j = 1:length(auxLFP(i).ch.Codes)
        foo = auxLFP(i).ch.Codes{j};
        trialRT(j) = foo(foo(:, 1) == 4, 2) - foo(foo(:, 1) == 62, 2);
    end
    
    % assign to data 
    sessRT{i} = trialRT(auxStruct.channelTrials{i});
    tmp = circshift(auxLFP(i).ch.Rewarded, 1); 
    rewarded{i} = tmp(auxStruct.channelTrials{i});
    chanInfo(i).sessRT    = sessRT{i}';
    chanInfo(i).rewarded  = rewarded{i};
end

auxStruct.sessRT   = sessRT;
auxStruct.rewarded = rewarded;

% ripple detector definitions 
auxStruct.skipRippleDetect = 1;
auxStruct.Fs               = 1000;
auxStruct.sdev             = 3;         % standard deviation cutoff
auxStruct.minDiff          = 40;        % min difference between individual peaks
auxStruct.lb               = 25;        % min ripple duration | vaz 2019
auxStruct.ub               = 200;       % max ripple duration
auxStruct.rippleBounds     = [80, 180]; % ripple band in nhps| logothetis 2012
auxStruct.smooth           = 0.015;     % smoothing factor
auxStruct.cutoffLength     = 15;        % min duration for >= cutoff
auxStruct.ripple_sd        = 6;

auxStruct.widebandBounds_supps = [2, 200];

% define filter for ripple fig
sm_inf = auxStruct.smooth*auxStruct.Fs;
sz    = sm_inf*auxStruct.ripple_sd;
x = linspace(-sz / 2, sz / 2, sz);             
gaussFilter = exp(-x .^ 2 / (2 * sm_inf ^ 2));  
gaussFilter = gaussFilter / sum (gaussFilter); 

auxStruct.gaussFilter = gaussFilter;

% neuron phase aligning-related definitions
nBins     = 11; 
radBounds = linspace(0, 2*pi, nBins); 

% correct for plot wrapping 
radBounds = mod(rad2deg(circshift(radBounds(1:(nBins) - 1) - pi, floor(nBins / 2))), 360);
radBounds(nBins) = 360;
auxStruct.radBounds = radBounds;
auxStruct.nBins = nBins;

% estimate bins using wavlength
wavLength = 249;
tBins     = [1:10:1500];
tBins     = tBins(tBins < 1400-wavLength);

auxStruct.tBins     = tBins;
auxStruct.grid_peak = 220:251;
auxStruct.neuron_theta_aligned_roi = 176:476;

% examines choice period up to 400 msec post choice; corresponds to timepoints
auxStruct.neuron_theta_averaging_roi = 1:1200;

% settings related to channel / unit info
[sessions, sessionsF, sessionsM, ~, ~, chanSessions, unitSessions] = extractSessions_v2(auxStruct.all_units_info);

auxStruct.sessions     = sessions;
auxStruct.sessionsF    = sessionsF;
auxStruct.sessionsM    = sessionsM;
auxStruct.chanSessions = chanSessions;
auxStruct.unitSessions = unitSessions;

auxStruct.rel_sr_tp = 21:141;     % timepoints of interest for the hex SR analysis
auxStruct.num_sessions_dat2 = 15; % num extracted sessions from dat2 

% peak / trough ROIs for vmpfc unit fig based on LFP values
auxStruct.vmpfc_roiStarts = [-266, -192, -127, -74]; % begin peak/trough within orig roi -345
auxStruct.vmpfc_roiEnds   = [-206, -132, -67, -14];  % end peak/trough within orig roi -315

auxStruct.fold_sym   = 6;
auxStruct.wavelength = 60;
auxStruct.cvals      = 3;
auxStruct.inputT     = auxStruct.startBin:auxStruct.endBin;
auxStruct.PE_roi     = auxStruct.inputT;

auxStruct.num_vmpfc_cells = 160;
auxStruct.nnmf_ripple_data_src = 'C:\Users\sebas\My Drive\Work\KennerleyLab\Code\Infogathering\cell_paper_revision\master\data\1st_stage_ripples\main\';
auxStruct.num_theta_timepoints = 476;
auxStruct.num_freqs = 14;

overlap_step                = 1;
smoothing_kernel            = 0.5; % sigma for kernel of 2d traces

auxStruct.overlap_step     = overlap_step;
auxStruct.smoothing_kernel = smoothing_kernel;
auxStruct.ripple_sd = 6;
auxStruct.maxIndx           = 900; % 
auxStruct.peaks             = [301:1200]; % corresponds to choice period within figures where max fr should be expected in the null 
auxStruct.ripple_fr_sup_roi = 150:200;
auxStruct.offset            = 300;

auxStruct.fr_baseline = 1:301;

auxStruct.nClusters  = 3;
auxStruct.SNR_THRESH = 5;
auxStruct.REPLICATES = 10;
auxStruct.logo_bounds = [22, 75, 80];

auxStruct.nnmf_ripple_fin_src = 'C:\Users\sebas\My Drive\Work\KennerleyLab\Code\Infogathering\cell_paper_revision\master\data\2nd_stage_ripples';
auxStruct.powerComp = 1; % power info for heatmap toggle

auxStruct.dataP = 'C:\Users\sebas\My Drive\Work\KennerleyLab\Code\Infogathering\cell_paper_revision\share_online\data\';
auxStruct.figP  = 'C:\Users\sebas\My Drive\Work\KennerleyLab\Figures\grids_project';

% params
auxStruct.alpha = 1;
auxStruct.gamma = 1;
auxStruct.theta = 10;
auxStruct.delta = 0;
auxStruct.zeta  = 0;


auxStruct.binom_thresh = .01;


end