
% Run heading to addpath and load networks
%{
open DesignTargets.m
%}

%% EA global, load data
result_file = 'subR1_EAgwtl2_pop2000.mat';
csv_Optim = 'outs-Nocut2_subR1_EAgwtl2_pop2000.csv';

% Go to Result_Check.m
%{
open Result_Check.m
%}

% Shape view
view(20,45);

%% run subdomain GA
% open GA_subspace.m

%% EA subdomain, load data

% subdomain EA 1
result_file = 'subR1_EAgwtl2_EAsub1disp_pop2000.mat';
csv_Optim = 'outs-Nocut2_subR1sub1disp_pop2000.csv';

% subdomain EA 2
result_file = 'subR1_EAgwtl2_EAsub2disp_pop2000.mat';
csv_Optim = 'outs-Nocut2_subR1sub2disp_pop2000.csv';

% % subdomain EA 3, no improvement
% result_file = 'subR1_EAgwtl2_EAsub3disp_pop2000.mat';
% csv_Optim = 'outs-Nocut2_subR1sub3disp_pop2000.csv';

% Go to Result_Check.m
%{
open Result_Check.m
%}

% Shape view
view(20,45);

