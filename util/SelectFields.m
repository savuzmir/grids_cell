function [allUnits] = SelectFields(units, fields, FileToLoad)

% This function selects individual fields of several units and combines them into structure that gets saved 
%INPUT:
%   units(CELL ARRAY):
%           Contains a cell array of all units that are to be loaded
%   fields (CELL ARRAY):
%           Contains a cell array of strings that has all the fields to be loaded in
%   AllNeurons (STRUCT ARRAY):
%           Contains a structure with information such as number of days and stimsets if it will be needed
%   FileToLoad (STRING):
%           Contains the file that is created in the end and warns that it already exists 
%OUTPUT:
%   out(MATRIX):
%           Matrix with indices of neurons to be sampled from the main matrix.
%           Each column represents indices for one region. 

% Sebastijan 26/3/19

if ~AlreadyExists(FileToLoad)  
    for i_unit = 1:length(units)
            tic
            currUnit = load(units{i_unit}, fields{:});
            allUnits(i_unit) = currUnit;
            sprintf('Unit %i needed %f', i_unit, toc)
    end
    
    
else
	fprintf('Already exists, loaded %s from folder.\n', FileToLoad)
    load(FileToLoad)
end
