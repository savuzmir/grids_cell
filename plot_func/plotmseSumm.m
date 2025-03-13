function h = plotmseSumm(varargin)

%  plotmse(varargin);
%
%  plots mean and standard error of ymat across dimension dim (default
%    1)


if nargin == 1
    ymat = varargin{1};
    colorinf = {'color', [84/255, 161/255, 140/255]};
elseif nargin == 2
    xsamp = varargin{1}';
    ymat = varargin{2};
    colorinf = {'color', [84/255, 161/255, 140/255]};
elseif nargin == 3
    xsamp = varargin{1}';
    ymat = varargin{2};
    colorinf = {'color', varargin{3}};
elseif nargin == 4
    xsamp = varargin{1}';
    ymat = varargin{2};
    colorinf = {'color', varargin{3}};
    fillAlpha = varargin{4};
end

if nargin == 1
    xsamp = 1:size(ymat, 1);
end

mn = ymat(:, 1)';
ub_se = ymat(:, 2)';
lb_se = ymat(:, 3)';


hold on;

% will fetch
% just color
filling = fill([xsamp fliplr(xsamp)], [ub_se fliplr(lb_se)], colorinf{2}, 'linestyle', 'none');
if nargin == 4
    alpha(filling, fillAlpha);
else
    alpha(filling, .2);
end

h=plot(xsamp, mn, 'LineWidth', 3, colorinf{:});
%plot(xsamp, ub_se, 'LineWidth',1, colorinf{:});
%plot(xsamp, lb_se, 'LineWidth',1, colorinf{:});

hold off;

end