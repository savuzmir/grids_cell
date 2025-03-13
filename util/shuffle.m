function [out, idx] = shuffle(inp)
% function [out, idx] = shuffle(inp)
% inp (VECTOR) [1xN or Nx1] of data to be shuffled
% out (VECTOR) of shuffled data with same dim
% ixd (VECTOR) of shuffle indices

% Generates a permutation of the data vector (inp)

idx = randperm(length(inp));
out = inp(idx);
end