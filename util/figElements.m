function figElements(figHandle, varargin)
% util function that makes prettier plots
% INPUT: title (STR), xlabel (STR), ylabel (STR), xlim (2VEC), ylim (2VEC), 
%        xticks (VEC), xticklabels (CELL), yticks (VEC), yticklabels (CELL), 
%        fontsize (DOUBLE), remove x ticks (INT), remove y ticks (INT)

    tit = varargin{1}; xAxis = varargin{2}; yAxis = varargin{3}; xLimit = varargin{4};
    yLimit = varargin{5}; xTickVals = varargin{6}; xTickLabels = varargin{7}; 
    yTickVals = varargin{8}; yTickLabels = varargin{9}; fontsize = varargin{10};
    removeXticks = varargin{11}; removeYticks = varargin{12};
  
  %  figure(figHandle);

    
    if all(~isempty(tit))  %removes assumptions individually
        title(tit);
    end
    
    if all(~isempty(xAxis))
    xlabel(xAxis);
    end
    
    if all(~isempty(yAxis))
    ylabel(yAxis);
    end
    
    if all(~isempty(xLimit))  
    xlim(xLimit);
    end
    
    if all(~isempty(yLimit))  
    ylim(yLimit);
    end
    
    if all(~isempty(xTickVals))  
    % add ticks, add tick info 
    xticks(xTickVals);
    end   
    
    if ~all(cellfun(@(xTickLabels) any(isempty(xTickLabels(:))), xTickLabels))
    xticklabels(xTickLabels);
    end
    
    if all(~isempty(yTickVals))  
    yticks(yTickVals);
    end
      
    if ~all(cellfun(@(yTickLabels) any(isempty(yTickLabels(:))), yTickLabels))
    yticklabels(yTickLabels);
    end
    
    if all(~isempty(fontsize))  
    set(gca, 'fontsize', fontsize);
    end
    
    if all(~isempty(removeXticks))  
    % add ticks, add tick info 
    set(gca,'xtick',[]);
    end   
    
	if all(~isempty(removeYticks))  
    % add ticks, add tick info 
    set(gca,'ytick',[])
    end   

end


