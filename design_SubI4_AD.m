
% Run heading to addpath and load networks
%{
open DesignTargets.m
%}

%% AD results load data
% Two types are slightly different on AD optimized pattern; 
% but with the same pattern, feeding (-1,1)---type2, or (-2mu,2-2mu)---type 3, 
% to the network gives very similar prediction.

% AD based on Type 2
result_file = 'subI4FlipX_AD0s.mat';
csv_Optim = 'outs-Nocut2_subI4FlipX_ADtype2_zeros.csv';

% Go to Result_Check.m
%{
open Result_Check.m
%}

% Shape view
view(-10,20);

%% run subdomain GA
% open GA_subspace.m

%% EA subdomain, load data

result_file = 'subI4FlipX_AD0s_EAsub1disp_pop2000.mat'; % RN11
csv_Optim = 'outs-Nocut2_subI4FlipX_AD0s_EARN11sub1disp_pop2000.csv';


% Go to Result_Check.m
%{
open Result_Check.m
%}

% Shape view
view(-10,20);

