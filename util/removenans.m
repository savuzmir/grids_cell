function out = removenans(inp, dim)
% removes nans in either rows or columns
% assumes, there is an equal number of them across rows/columns

if dim == 1
    inp = inp(~isnan(inp(:, 1)), :);
else
    inp = inp(:, ~isnan(inp(1, :)));
end
out = inp;
end