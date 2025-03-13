function [out] = demean(mainIn, demeanIn, demeanType, auxStruct, field, fieldPre, onlyBin)
tic
Units = 1:length(mainIn);

% 3 types of demeaning/avging
% 3. demean + avg within a specific time window 


% Preprocess the matrix you want according to specifications 
windowSize =    auxStruct.windowSize; 
slideWidth =    auxStruct.slideWidth;
epochDur =      auxStruct.epochDur;

if ~onlyBin
    if demeanType == 1 | demeanType == 2 | demeanType == 3
    preEvent =      auxStruct.preEvent;
    slSizeStart =   round(auxStruct.slSizeStart/slideWidth);
    slSizeEnd =     round(auxStruct.slSizeEnd/slideWidth);
    preEventBins =  round((preEvent/slideWidth)) + 1;
    end
end

% perform binning and smoothing, preserve structure type 

totBins=        (epochDur/slideWidth) - ((windowSize/slideWidth)-1);
bins =          [1:slideWidth:(epochDur - slideWidth + 1)];

smoothedData = struct;

if onlyBin

    for u = Units                        % just bins your data according to predefined criteria 
        for currBin= 1:floor(totBins)
            smoothedData(u).Unit(:, currBin)=nanmean(mainIn(u).(field{1})(:, bins(currBin):(bins(currBin) + (windowSize - 1))), 2);
        end
    end

else
    
    for u = Units                       
        for currBin= 1:floor(totBins)
            smoothedData(u).Unit(:, currBin)=nanmean(mainIn(u).(field{1})(:, bins(currBin):(bins(currBin) + (windowSize - 1))), 2);
        end
    end
    
    switch demeanType 
        
        case 1           % just demean according to a preEvent period                                                                                                                                      

        for u = Units
            FR = smoothedData(u).Unit;
            preEventMean =            nanmean(FR(:, 1:preEventBins), 2);                                                                                                   
            LockedPostEventMean =     FR - preEventMean;                                    
            smoothedData(u).Unit = LockedPostEventMean;
        end
        
        case 2           % just demean according to a slice within the preEvent period                                                                                                                                      

        for u = Units
            FR = smoothedData(u).Unit;
            preEventMean =            nanmean(FR(:, slSizeStart:slSizeEnd), 2);                                                                                                   
            LockedPostEventMean =     FR - preEventMean;                                    
            smoothedData(u).Unit = LockedPostEventMean;
        end

        case 3          % demean according to preEvent period + take the average FR of a slice 

        %  first smooth, then take a bin and average all firing within that bin 
        sliceStart =    preEventBins + slSizeStart;                                              % We select a timewindow between slSizeStart and slSizeEnd msec post event
        sliceEnd =      preEventBins + slSizeEnd;                                                
 
        for u = Units
            FR = smoothedData(u).Unit;
            preEventMean =            nanmean(FR(:, 1:preEventBins), 2);                                                                                                   
            LockedPostEventMean =     FR - preEventMean;                                    
            postEventMean =           nanmean(LockedPostEventMean(:, sliceStart:sliceEnd), 2);
            smoothedData(u).Unit = postEventMean;

        end    

        case 4          % demean according to another struct (e.g. fixation period)
        windowSizePre =    auxStruct.windowSizePre; 
        slideWidthPre =    auxStruct.slideWidthPre;
        epochDurPre =      auxStruct.epochDurPre;

        totBins=        (epochDurPre/slideWidthPre) - ((windowSizePre/slideWidthPre)-1);
        bins =          [1:slideWidthPre:(epochDurPre - slideWidthPre + 1)];

        smoothedDataPre = struct;

        for u = Units
            for currBin = 1:floor(totBins)
                smoothedDataPre(u).Unit(:, currBin)=nanmean(demeanIn(u).(fieldPre{1})(:, bins(currBin):(bins(currBin) + (windowSizePre - 1))), 2);
            end
        end


        for u = Units
            FR = smoothedData(u).Unit;
            preFR = nanmean(smoothedDataPre(u).Unit, 2);           
            LockedPostEventMean =     FR - preFR;                                    

            smoothedData(u).Unit = LockedPostEventMean;
        end
        
        
        case 5                   % demean according to a slice of another struct (e.g. fixation period)
            
            windowSizePre =    auxStruct.windowSizePre;
            slideWidthPre =    auxStruct.slideWidthPre;
            epochDurPre =      auxStruct.epochDurPre;
            
            totBins=        (epochDurPre/slideWidthPre) - ((windowSizePre/slideWidthPre)-1);
            bins =          [1:slideWidthPre:(epochDurPre - slideWidthPre + 1)];
            
            smoothedDataPre = struct;
            
            sliceStart     =    round(auxStruct.slSizeStartPre/slideWidthPre);
            sliceEnd       =    round(auxStruct.slSizeEndPre/slideWidthPre);
               
            for u = Units
                for currBin = 1:floor(totBins)
                    smoothedDataPre(u).Unit(:, currBin)=nanmean(demeanIn(u).(fieldPre{1})(:, bins(currBin):(bins(currBin) + (windowSizePre - 1))), 2);
                end
            end
            
            for u = Units
                FR = smoothedData(u).Unit;
                preFR = nanmean(smoothedDataPre(u).Unit(:, sliceStart:sliceEnd), 2);
                [nRow, nCol] = size(FR);
                LockedPostEventMean =     FR - repmat(preFR, [1, nCol]);
                
                smoothedData(u).Unit = LockedPostEventMean;
            end
            
           
    end
end

out = smoothedData;
toc
end


