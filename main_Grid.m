clc;clear;

%% preliminaries

% set base path
addpath(genpath('folder with main func'));

%% run fig functions
% individual functions load relevant files to create figures as needed
% some files are large so could cause ram 'out of memory' issues depending on pc specs

% generate figure 1 - behavior + value
figure1_BehaviourValue();
figure1_supplements_warping();

% generate figure 2 - grid code
figure2_lfpGrid();
figure2_lfpGrid_supplements();

% generate figure 3 - grid orientation and consistency 
figure3_gridRealignment();

% generate figure 4 - cell-related results
figure4_GridNeurons();
figure4_GridNeurons_supplements(); 
ofc_post_revision_supplements();

% generate figure 5 - ripple related results
figure5_ripples();
figure5_rippleSupplements(); 
nnmf_example_supplements(); 




