function [ranks]=calc_rank(inp)

ranks = nan(size(inp, 1), 1);

uniqueRanks = unique(inp);
rmIndx = find(uniqueRanks~=0);
uniqueRanks = uniqueRanks(rmIndx);

rankedVals = [];

for i = 1:size(uniqueRanks)
    insertIndx = find(uniqueRanks(i) == inp);
    ranks(insertIndx) = i;
end

end