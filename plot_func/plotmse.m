function h = plotmse(varargin)
% using this function

% ==========================
% 1. basic 
% use:
% smooth_bool   = 0;
% smooth_length = 0;
% data          = randn(50, 1500);
% plotmse(data, [smooth_bool, smooth_length]);

% explained: 
% plotmse(data, [smooth_bool, smooth_length])
% data:          [nFeatures x nTimepoints]
% smooth_bool:   [0/1]. Do you want to smooth along nTimepoints?
% smooth_length: [int {1, nTimepoints-1}]. How many timepoints along nTimepoints to smooth over?

% ==========================
% 2. coloring the trace
% use:
% color         = [0.5, 0.5, 0.5];
% plotmse(data, color, [smooth_bool, smooth_length])

% explained:
% plotmse(data, color, [smooth_bool, smooth_length])
% color:         [matlab string OR RGB vector]. Colors the mean and SEM. 

% ==========================
% 3. modifying transparency of SEM
% use:
% sem_alpha     = 0.5;
% plotmse(data, color, [smooth_bool, smooth_length], sem_alpha) 

% explained:
% plotmse(data, color, [smooth_bool, smooth_length], sem_alpha)
% sem_alpha:     [float {0, 1}]. Controls data alpha 

% ==========================
% 4. modifying everything
% use:
% xDim          = -300:1:1199;
% color2        = 'r';
% lineWidth     = 0.5;
% plotmse(xDim, data, color1, color2, [smooth_bool, smooth_length], sem_alpha, lineWidth)

% explained:
% xDim: [1 x nTimepoints]. Customizes X axis.
% color1:       [matlab string OR RGB vector]. Color of the mean.
% color2:       [matlab string OR RGB vector]. Color of the SEM.
% lineWidth:    [float {0.001, inf-1}]. Thickness of mean.

% ==========================
% 5. using the legend with several averages
% use:
% colMat       = {'r', 'g', 'b'};
% leg_arr      = [];
% data         = randn(50, 1500, 3);
% for i = 1:3
% leg_arr(i) = plotmse(data(:, :, i), colMat{i}, [0, 1]);
% end
% legend(leg_arr(:), {'a', 'b', 'c'})

% random cool colors
% cols = {[230/255, 159/255, 1/255],   [86/255, 180/255, 233/255], ...
%         [1/255,   158/255, 115/255], [240/255, 228/255, 66/255], ...
%         [1/255,   114/255, 178/255], [213/255, 94/255, 1/255], ...
%         [204/255, 121/255, 167/255]};


cols =  {[148/255, 149/255, 153/255], ...
         [148/255, 149/255, 153/255], ...
         [148/255, 149/255, 153/255], ...
         [148/255, 149/255, 153/255], ...
         [84/255, 161/255, 140/255]};

if nargin == 2
    ymat = varargin{1};
    colorinf = {'color', [0.3294, 0.6314, 0.5490]};
	smoothing = varargin{2};
elseif nargin == 3
    ymat = varargin{1};
    colorinf = {'color', varargin{2}};
    smoothing = varargin{3};
elseif nargin == 4
    ymat = varargin{1};
    colorinf = {'color', varargin{2}};
    smoothing = varargin{3};
    fillAlpha = varargin{4};
elseif nargin >= 5
    xDim = varargin{1};
    ymat = varargin{2};
    colorinfM = {'color', varargin{3}, 'linewidth', varargin{7}};
    colorinfS = {'color', varargin{4}, 'markeredgecolor', 'none'};
    smoothing = varargin{5};
    fillAlpha = varargin{6};    

end

xsamp = 1:size(ymat, 2);

higherDim = size(ymat, 3);

if higherDim > 1 && ~smoothing(1)
    mn = squeeze(nanmean(ymat, 1))';
    se = squeeze(nanstd(ymat, [], 1)./sqrt(size(ymat,1)))';
elseif higherDim == 1 && ~smoothing(1) 
    mn = nanmean(ymat, 1);
    se = nanstd(ymat, [], 1)./sqrt(size(ymat(find(~isnan(ymat(:, 1))), :), 1));
end

if smoothing(1)
    
    if higherDim > 1
        for cdim = 1:higherDim
            mn(cdim,:) = smooth(nanmean(ymat(:, :, cdim), 1), smoothing(2))';
            se(cdim,:) = smooth(nanstd(ymat(:, :, cdim), [], 1)./sqrt(size(ymat(find(~isnan(ymat(:, 1))), :), 1)), smoothing(2))';
        end
    else
        
        mn = smooth(nanmean(ymat, 1), smoothing(2))';
        se = smooth(nanstd(ymat, [], 1)./sqrt(size(ymat(find(~isnan(ymat(:, 1))), :), 1)), smoothing(2))';
        
    end
end

ub_se = mn + se;
lb_se = mn - se;
hold on;

% will fetch
% just color
if higherDim > 1
    for cdim = 1:higherDim
        colorinfS = {'color', cols{cdim}, 'markeredgecolor', 'none'};
        
        if nargin >= 5
            filling = fill([xDim fliplr(xDim)], [ub_se(cdim, :) fliplr(lb_se(cdim, :))], colorinfS{2}, 'linestyle', 'none');
        else
            filling = fill([xsamp fliplr(xsamp)], [ub_se(cdim, :) fliplr(lb_se(cdim, :))], colorinf{2}, 'linestyle', 'none');
        end

        if nargin >= 4
            alpha(filling, fillAlpha);
        else
            alpha(filling, .2);
        end
    end
    
else
    if nargin >= 5
        filling = fill([xDim fliplr(xDim)], [ub_se fliplr(lb_se)], colorinfS{2}, 'linestyle', 'none');
    else
        filling = fill([xsamp fliplr(xsamp)], [ub_se fliplr(lb_se)], colorinf{2}, 'linestyle', 'none');
    end
    
    if nargin >= 4
        alpha(filling, fillAlpha);
    else
        alpha(filling, .2);
    end
    
end

% plot means
if higherDim > 1
    for cdim = 1:higherDim
        
        colorinfM = {'color', cols{cdim}, 'linewidth', varargin{7}};
        
        if nargin >= 5
           h = plot(xDim, mn(cdim, :), colorinfM{:});
        else
           h = plot(xsamp, mn(cdim, :), 'LineWidth',2, colorinf{:});
        end
    end
    
else
    if nargin >= 5
        h = plot(xDim, mn, colorinfM{:});
    else
        h = plot(xsamp, mn, 'LineWidth',2, colorinf{:});
    end
end

end