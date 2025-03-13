function out = normalize_bound(inp, lb, ub, varargin)
% function out = normalize_bound(inp, lb, ub, varargin)
% if it's a matrix we assume a trials x time format

mat_size = size(inp);

% if its a matrix, we normalize within trial - assumes trials x time format
if mat_size(1) > 1 & mat_size(2) > 1

    ma = max(inp, [], 2);
    mi = min(inp, [], 2);

    out = (ub - lb) *((inp - mi) ./ (ma - mi)) + lb;

else
    if isempty(varargin)
        ma = max(inp(:));
        mi = min(inp(:));
    else
        ma = varargin{1};
        mi = varargin{2};
    end
    
    out = (ub - lb) *((inp - mi) ./ (ma - mi)) + lb;
end

end