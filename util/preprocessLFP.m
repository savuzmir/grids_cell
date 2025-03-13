function chanInfo = preprocessLFP(data, powerTransform, freqSettings, Fs, rel_voices) 
% This function notch filters (removing 50 and harmonics + interpolates them with neighbouring frequencies) and bandpassses data between 0.1 and 200 Hz
% if powerTransform is set to 1, it also performs a morlet wavelet convolution

% inputs:
%   data (MATRIX: ntrials x ntime) for a channel  
%   powerTransform (BOOL) - do you want to perform the morlet transform 
%   freqSettings ([1x2]) - min and max bandpass width
%   Fs - sampling freq
%   rel_voices - granularity of power computation

% outputs:
%       chanInfo (STRUCT) with fields 
%           lfp (ntrials x ntime) notchfiltered and bandpassed lfp signal
%           lfpPower (ntrials x ntime x frequencies) power convolved with a morlet wavelet
                                                                        
% parameters
% bandpass we are using for the data preprocessing 
bandpassWidth = freqSettings;
lineRemoval = [50, 100, 150; 1 1 1; 3 3 3];                                 

% voices per octave - how dense is the sampling of frequencies
voices = rel_voices; % 20; 
% frequency limits we are interested in
freqLims = freqSettings;
                                           
% extract lfp signal
currChannel = data;

% time length 
samples     = size(currChannel, 2);


% generate filter settings for frequency extraction
fb = cwtfilterbank('SamplingFrequency', Fs, 'SignalLength', samples, 'FrequencyLimits', freqLims, 'wavelet', 'amor', 'VoicesPerOctave', voices);

% generate container
chanInfo              = struct;
chanInfo.header       = {'Trials x Time x Frequencies'};
chanInfo.lfp          = [];
chanInfo.frequencies  = [];
chanInfo.coi          = [];
chanInfo.lfpMagnitude = [];
chanInfo.settings     = fb;

% number of trials within session 
trLen = size(currChannel, 1);

for tr = 1:trLen
    
    tmp = currChannel(tr, :);
    
    tmp = bandpassfilter(ft_preproc_dftfilter(tmp, Fs, lineRemoval(1, :), 'Flreplace', 'neighbour', 'Flwidth', lineRemoval(2, :), 'NeighWidth', lineRemoval(3, :)), Fs, bandpassWidth);
    chanInfo.lfp(tr, :) = tmp;
    
end

% if power transform
if powerTransform
    for tr = 1:trLen

        tmp = chanInfo.lfp(tr, :);
        [vals, freq, coi] = cwt(tmp, 'FilterBank', fb); % morlet wavelet of our time window of interest  
        chanInfo.lfpMagnitude(tr, :, :) =  abs(vals)';
    end
    % needs only one writing per channel
    chanInfo.frequencies = freq;
    chanInfo.coi = coi;
end

end





