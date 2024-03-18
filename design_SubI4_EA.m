
% Run heading to addpath and load networks
%{
open DesignTargets.m
%}

%% Check global EA resuls

% % EAg, non-weighted loss
% result_file = 'subI4FlipX_EAg_pop2000.mat';
% csv_Optim = 'outs-Nocut2_subI4FlipX_EAconvRN22_pop2000.csv';

% EAg, weighted loss 2
result_file = 'subI4FlipX_EAgwtl2_pop2000.mat';
csv_Optim = 'outs-Nocut2_subI4FlipX_EAconvRN22_pop2000_wtloss2.csv';

% Go to Result_Check.m
%{
open Result_Check.m
%}

% Shape view
view(-10,20);

%% run subdomain GA
% open GA_subspace.m

%% Check subdomain EA results

% % subdomain EA 1
% load subI4FlipX_EAg_EAsub1disp_pop2000.mat
% csv_gaOptim = 'outs-Nocut2_subI4FlipX_EAg_EAsub1disp_pop2000.csv';

% subdomain EA 1
result_file = 'subI4FlipX_EAgwtl2_EAsub1disp_pop2000.mat';
csv_Optim = 'outs-Nocut2_subI4FlipX_EAgwtl2_EAsub1disp_pop2000.csv';

% subdomain EA 2
result_file = 'subI4FlipX_EAgwtl2_EAsub2disp_pop2000.mat';
csv_Optim = 'outs-Nocut2_subI4FlipX_EAgwtl2_EAsub2disp_pop2000.csv';

% Go to Result_Check.m
%{
open Result_Check.m
%}

% Shape view
view(-10,20);

