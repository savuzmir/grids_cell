function out = nanlength(inp)

inp = length(inp(~isnan(inp)));
out = inp;

end