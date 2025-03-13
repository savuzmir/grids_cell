function out = find_rows(inp1, inp2, repeated)
% function out = find_rows(inp1, inp2, repeated)
% find rows of INP1 (COLUMN VECTOR) with values of INP2 (COLUMN VECTOR) 
% OUT (VECTOR) are the rows of INP1 with containing values of INP2

% repeated can be used if you expect several instances of a given element in inp2 to be present in inp1 -
% out is then a cell for each instance of inp2


if strcmp(repeated, 'yes')
    out = {};
else
    out = [];
end

if size(inp2, 1) == 1 && size(inp2, 2) ~= 1
    inp2 = inp2';
end

if size(inp1, 1) == 1 && size(inp1, 2) ~= 1
    inp1 = inp1';
end

for i = 1:size(inp2, 1)
    
    try % if this doesn't work, then there are several 
        if strcmp(repeated, 'yes')
            out{i} = find(inp1 == inp2(i));
        else
            out(i) = find(inp1 == inp2(i));
        end
    catch
        out = out + (inp1 == inp2(i))';
    end
   
end

if strcmp(repeated, 'yes')
    out = cell2mat(out');

end
