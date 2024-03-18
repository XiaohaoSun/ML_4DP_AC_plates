
% Run heading to addpath and load networks
%{
open DesignTargets.m
%}

global outIntui

%% GD, load data
% Two types are slightly different on AD optimized pattern; 
% but with the same pattern, feeding (-1,1)---type2, or (-2mu,2-2mu)---type 3, 
% to the network gives very similar prediction.

% AD based on Type 2
result_file = 'AD_SubR2_ones.mat';
csv_Optim = 'outs-Nocut2_subR2_ADtype2_ones.csv';

result_file = 'AD_SubR2_zeros.mat';
csv_Optim = 'outs-Nocut2_subR2_ADtype2_zeros.csv';

result_file = 'AD_SubR2_neutral.mat';
csv_Optim = 'outs-Nocut2_subR2_ADtype3_neutral.csv';

result_file = 'AD_SubR2_random.mat';
csv_Optim = 'outs-Nocut2_subR2_ADtype3_random.csv';

% Go to Result_Check.m
%{
open Result_Check.m
%}

% Shape view
view(0,35);

%% EA global, load data
result_file = 'subR2_EAgwtl2_pop2000.mat';
csv_Optim = 'outs-Nocut2_subR2_EAgwtl2_pop2000.csv';

% Go to Result_Check.m
%{
open Result_Check.m
%}

% Shape view
view(0,35);

