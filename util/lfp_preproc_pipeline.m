%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this pipeline is designed such that everything is computed from the raw datafiles 
% as intermediary steps, containers with raw or preprocessed LFP and auxiliary information can be stored
% note that pathing / folders were specific to original pc on which data was processed and have been changed so may not work 

clc; clear;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% join individual LFP datasets into one file(ContainerLFP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                           
baseDir  = ''; 
rawFiles = [baseDir, '\preprocessed_units\'];                                                              % folder with preprocessed units

folder0  =  [baseDir, '\0_raw_lfp\'];
folder0b =  [baseDir, '\0_raw_lfp_reduced\'];
folder0c =  [baseDir, '\0_raw_FR\'];
folder1  =  [baseDir, '\1_bandpass_separated\'];                                                            % is used at various stages for saving 
folder2  =  [baseDir, '\2_morlet_transform\']; 
folder3  =  [baseDir, '\3_binned\'];                                  

codepath = '';

addpath(codePath)
addpath([codePath, '\util'])                                                                                % add all functions to current path

cd(rawFiles)                                        

allFiles = dir();                                                                                           % store the preprocessed units 
allFiles = {allFiles.name};                                                                         
allFiles = {allFiles{3:end}};                                                                               % the linux distro finds empty info that needs to be removed 

% File to load and main fields can have different names, hence the representation as two different variables. 

% for FR 
FileToLoad = {'ResponseMatrix'}; % if file already exists, it prevents it from running the function
mainFields = {'ResponseMatrix'};

% Save preprocessed individual units into ContainerLFP for individual task-relevant matrices
for file = 1:length(FileToLoad)
    
    AuxiliaryLFP = struct;                                                                                     % containerLFP has main LFP data
                                                                                                               % auxiliaryLFP has additional, task-relevant data
    
    ContainerFR = SelectFields(allFiles, mainFields(file), FileToLoad{file});                                  % SelectFields loops through raw files and adds lfp matrices from mainFields into ContainerLFP
    
    % independently save auxiliary information - using the same principle as above
    % simultaneous trials
    auxiliaryFields = {'Prob_top', 'Left_pay', 'Left_prob', 'Right_pay', 'Right_prob', 'spatial_values', 'rnkSpatVal', ...
                      'sdev', 'men', 'ChEVL', 'ChEVR', 'ChML', 'ChMR', 'ChPL', ...
                      'ChPR', 'UnEVL', 'UnEVR', 'UnML', 'UnMR', 'UnPL', 'UnPR', ...
                      'Rewarded', 'brain_region', 'sb', 'Codes', 'plexCodes', ...
                      'postCue', 'postRew', 'postFix', 'postResp', ...
                      'preCue', 'preFix', 'preResp', 'preRew'};
    
    if file == 2        % Save it only on the first run if you are saving multiple matrices in one for loop, should be 1 
        
        % infogathering
        auxiliaryFields = {'first_picture', 'second_picture', 'third_picture', 'Right_EV', 'Left_EV', 'brain_region', 'brainer', ...
            'sdev', 'men', 'Codes', 'condition', 'horz', 'Left_pay', 'Left_prob', ...
            'LM_seen', 'LP_seen', 'P_first', 'P_second', 'P_third', 'RM_seen', 'RP_seen', ...
            'spatial_values', 'vert', 'sb', 'Plex_codes', 'Left_choice'...
            'post_cue', 'postFix', 'postResp', 'pre_cue', ...
            'preFix', 'preResp'};
        
        
        auxFile = 'AuxiliaryData.mat';
        AuxiliaryLFP = SelectFields(allFiles, auxiliaryFields, auxFile);
        
        
        % At this point all the relevant files can be saved 
        saveStep0 = 1;

        if saveStep0
            full = [folder0, mainFields{file}];
            save(full, 'ContainerLFP', 'AuxiliaryLFP', '-v7.3');
        end

    else
        
        if saveStep0
            full = [folder0c, mainFields{file}];
            save(full, 'ContainerFR', '-v7.3');
        end
        
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% reduce repeated channels into a smaller structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('all_units_info.mat')

[correctCells, sessions] = findUniqueChannels(allFiles, all_units);             % Retrieves unique channel identifiers from all units
Units = 1:size(correctCells, 1);                                                % correctCells corresponds to channels with corresponding indices to units 
                                                                                % Units is a generic loop counter with 1:nr channels

cd('\0_raw_lfp');
allFiles = dir();
allFiles = {allFiles.name};
allFiles = {allFiles{3:end}};

saveRegressor = 0;                                                              % Do you want to save the angle regressor?
                                                                                
MatrixFields            = {'ResponseMatrixLFP'};                                                                         
                                                                                
for file = 1:length(MatrixFields)
    
    load(allFiles{file});   
    % we select the original containers obtained above and reduce them correctly

    currLFP                 = struct;
    currLFP(1).ch           = [];
    auxLFP                  = struct;
    
    for u = Units
         currLFP(u).ch       = ContainerLFP(correctCells(u)).(MatrixFields{file});     % from the original container with all units, we only select unique channels
        
        if file == 1                                                                 
            auxLFP(u).ch        = AuxiliaryLFP(correctCells(u));                      % The same is done for AuxiliaryLFP
        end
    end
   
    if saveRegressor
        %==================================================================
        %%%============================================= generate regressor
        %==================================================================
        RegressorContainer                 = generateAngleRegressor(AuxiliaryLFP);              % from the full AuxiliaryLFP we can generate the qf angle regressor for our data
        RegressorContainer(1).correctCells = correctCells;                                      % out of convenience we save the correctCells for later use
    end
    
    % At this point all the relevant files can be saved
    saveStep1 = 1;
    
    if file == 1
        if saveStep1
            full = [folder0b, MatrixFields{file}];                                                         
            save(full, 'currLFP',      '-v7.3');
            save('auxLFP');
            save('correctCells');
            save('sessions');
        end
        
    else
        full = [folder0b, MatrixFields{file}];                                                         
        save(full, 'currLFP', '-v7.3');
    end
    
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bandpass and notch filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd([baseDir, '\0_raw_lfp_reduced']);
allFiles = dir();
allFiles = {allFiles.name};
allFiles = {allFiles{3:end}};

yes_nans = 0;
Fs = 1000;                                      % sampling frequency
bands = {'highThetaB', 'betaB', 'lowGammaB', 'highGammaB'};
bandFreqs = [5 8; 12 30; 40 60; 61 100];       % these values correspond to Maidenbaum et al. (2018) PNAS

for file = 1:length(allFiles)
    
    load(allFiles{file});
  
    % At this point, the cleaned individual bands can be saved as necessary
    saveStep2 = 1;
    
    currFile = {'FixMatrixLFP', 'ResponseMatrixLFP', 'RewardMatrixLFP'};

    Units = 1:size(currLFP, 2);
    
    for u = Units
        trialLen(u) = size(currLFP(u).ch, 1);
    end
    
    timeLen = size(currLFP(u).ch, 2);
    
    for b = 1:length(bands)
        tmp = cell(1, 1, length(Units));
        
        % define relevant bands | allFreq will save the data structure with all bands combined from 1 to 200 Hz
        currLFP(1).(bands{b})    = [];
        
        parfor u = Units
            tic
            
            % only necessary if there are nans in the data
            if yes_nans
                
                tmpTrial = nan(trialLen(u), timeLen);
                
                for tr = 1:trialLen(u)
                    try
                        tmpTrial(tr, :) = bandpassfilter(ft_preproc_dftfilter(currLFP(u).ch(tr, :),  Fs, [50, 100, 150, 200]), Fs, bandFreqs(b, :));  % Fieldtrip's bandpass filter (enter: electrode x time)
                    catch
                    end
                end
                tmp{:, :, u} = tmpTrial;
            else
                
                tmp{:, :, u} = bandpassfilter(ft_preproc_dftfilter(currLFP(u).ch,  Fs, [50, 100, 150, 200]), Fs, bandFreqs(b, :));
                
            end
            toc
        end
        
        for u = Units
            currLFP(u).(bands{b}) = tmp{:, :, u};
        end
        
        if saveStep2
            full = [folder1, currFile{file}, '_', bands{b}];                                                        
            save(full, 'currLFP', '-v7.3');
        end
        currLFP = rmfield(currLFP, bands{b});
        
    end
    
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% now we perform a morlet  transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load in relevant auxiliary datasets from before
cd('1_raw_separated');          
load('RegressorContainer.mat');                                           

cd(baseDir)
load('auxLFP.mat');

load_from_file = 1;                                                        % if in the previous step you saved the file and want to load it
% then let it be 1, otherwise set to 0

if load_from_file
    cd('\1_bandpass_separated');                                % this loads in raw bandpassed and notch-filtered LFP saved above
    allFiles = dir();                                           % if from the previous steps there were files that were saved
    allFiles = {allFiles.name};                                 % in the 1_raw_separated, the following sequence can be run as is
    allFiles = allFiles(3:end);
else
    allFiles = 1;                                                          % set this to the number of fields you want to process
end

saveStep3 = 1;                                                    % if you to save a power transform before continuing

% name can be programatically changed by coding it as a variable
fileCtr      = 1;
onlyBin      = 1;                                                          % we only want to bin information in the demean function, as opposed to subtracting values w.r.t some preEvent period

inField      = {'MThighThetaB', 'MTbetaB', 'MTlowGammaB', 'MThighGammaB'};       % denotes the fields we are going to index. Done due to same reasoning as Variable
currFile = {'ResponseMatrixLFP'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% these values are specific to the sequence you run, the order needs to be based on the infield
indexVal    =   '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for file = 1:length(allFiles)                                           
    
    currBand = inField(indexVal(fileCtr));                              % tells us which band we are currently indexing.
    
    if load_from_file
        load(allFiles{file});                                           % file to be loaded in
        disp(sprintf('Successfully loaded %s', allFiles{file}))
    end
            
    inField      = {'MThighThetaB', 'MTbetaB', 'MTlowGammaB', 'MThighGammaB'};       % denotes the fields we are going to index. Done due to same reasoning as Variable
    folder2 = [baseDir '\2_morlet_transform\'];
    prefix =  'MT_';

    currLFP = powerTransform(currLFP, indexVal(fileCtr), 1);
    
    if saveStep3
        
        full = [folder2, prefix, allFiles{file}];   % has some hardcoded structs for magnitude, probability, chosen unchosen, etc.
        tmp = currLFP;
        currLFP = rmfield(currLFP, {bands{indexVal(fileCtr)}});
        currLFP
        disp(sprintf('About to save %s', allFiles{file}))
        save(full, 'currLFP', '-v7.3');
        disp('Successfully saved')
        currLFP = tmp;                      % we take back the original and delete just the last transform but not the band
        currLFP = rmfield(currLFP, {inField{indexVal(fileCtr)}});
        disp('Successfully cleared and ready for next iteration')
    end

    fileCtr = fileCtr + 1;
    disp(sprintf('fileCtr increment: %i', fileCtr))
end

%==================================================================
%%%=================================================== bin the data
%==================================================================

folders = {'\2_morlet_transform'};

startSearch =   [100];                                                  
endSearch   =   [950];                                 
fileOnset   =   [650]; % centers the artefacual trial search                                                                                                                   
                                                                       
% the output of this section becomes the final Y and X variables

auxStruct.onset         = fileOnset(file);
auxStruct.startSearch   = startSearch(file);                             % a timewindow used to exclude trials
auxStruct.endSearch     = endSearch(file);                               % with more than 10% timepoints which have a value of
                                                                         % 3 SD above or below the mean at the corresponding timepoint (chosen arbitrarily)
                                                                         % exludes from 0 to 4% percent of trials in all channels
                                                                         % Preprocess the matrix you want according to specifications
auxStruct.windowSize    = 100;                                           % what is the window size of interest for binning
auxStruct.slideWidth    = 10;                                            % what is the sliding window used for binning

saveStep4         = 1;
removeNoisyTrials = 0; % need to be careful at what step they are removed to ensure trials line up 
inField           = {'MThighThetaB', 'MTbetaB', 'MTlowGammaB', 'MThighGammaB'};  

for type = 1:length(folders)
    
    
    if load_from_file
        cd([baseDir, folders{type}]);                                       % this loads the morlet transform of all possible datasets
        allFiles        =  dir();
        allFiles        = {allFiles.name};
        allFiles        = allFiles(3:end);
    else
        allFiles        = 1;                                                % set this to the number of fields you want to process
    end
    
    VariableY	 = {'Y'};                                                   % name of our Y variable
    
    
    for file = 1:length(allFiles)                                           % This part seems to run into issues if it's being parfor-ed
        
        currBand = inField(indexVal(fileCtr));                              % tells us which band we are currently indexing.
        
        if load_from_file
            load(allFiles{file});                                           % file to be loaded in
            disp(sprintf('Successfully loaded %s', allFiles{file}))
            auxStruct.epochDur   =    size(currLFP(1).ch, 2);               % always corresponds to epoch duration
        end
       
        yBinned = demean(currLFP, [], [], auxStruct, currBand, [], onlyBin);
        
        for u = Units
            currLFP(u).(VariableY{:})   = yBinned(u).Unit;
            currLFP(u).aux              = auxLFP(u).ch;
        end

        [~, NoisyTrials] = getMax_Lag(currLFP, auxStruct, 'ch', 1);
        
        if removeNoisyTrials
            for u = Units
                trLen = (1:size(currLFP(u).(VariableY{1}), 1))';
                currLFP(u).(VariableY{:})   = currLFP(u).(VariableY{:})(setdiff(trLen, NoisyTrials(u).trialNums), :);
                currLFP(u).removeTrials     = NoisyTrials(u).trialNums';
                currLFP(u).aux              = auxLFP(u).ch;
            end
        end
        
        currLFP = rmfield(currLFP, {currBand{1}, 'ch'});
        
        if saveStep4
            full = [folder3, int2str(auxStruct.windowSize), 'Msec_', allFiles{file}];
            save(full, 'currLFP', '-v7.3');
        end
     fileCtr = fileCtr + 1;
    end
end

%% run hex analysis on these files



