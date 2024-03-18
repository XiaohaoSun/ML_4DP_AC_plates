
% result_file, stored by the user.
% Some existing results are stored in "mat_designs".
result_path = fullfile('mat_designs', result_file);

%% Read results

% Results from GD
load(result_path,'dlin','dlnetxy','dlnetz');
ga_x_mod = reshape(round_of_dlin(dlin)*0.5+0.5, 1,[]);

% Results from EA global
load(result_path,'ga_x');
ga_x_mod = ga_x;

% Results from EA subdomain
load(result_path,'ga_x','subspace','subspaceInd','gaSub_x','gaSub_pops','gaSub_out');
ga_x_mod = ga_x; ga_x_mod(subspaceInd) = gaSub_x;

%% FE verification
% =========================================================================
% If FE verification is needed, generate ga_insOptim.
ga_insOptim = reshape(ga_x_mod,[15,15,2]); % shape 15x15x2
ga_insOptim = [ga_insOptim(:,:,1); ga_insOptim(:,:,2)]; % shape 30x15

% "ga_insOptim" is copied to the voxel_input.txt,
% then go there to run "Exp2D_nocut_BC3_cae.py" for FE simulation.
% This will generate a csv file. Give its name to "csv_Optim" below:
csv_Optim = '';
% =========================================================================

%% Check patterns

% show subdomain being optimized
visualizeSubdomain(subspace);

% show "insIntui" if it exists
visualizeDesign(insIntui);

% show optimized pattern
visualizeDesign(ga_x_mod);

%% Check shapes

paras_nets = {netxy,netz, mu,mux,muy,muz,sig,sigx,sigy,sigz};

% Evaluation based on point pair distance
% ShapeDistP2P('GA','convBC', outIntui,ga_x_mod,csv_Optim, paras_nets{:})
% errorTable('GA','convBC', outIntui,ga_x,csv_Optim,insIntui, paras_nets{:})

% Evaluation based on point-to-surface distance
ShapeDist('GA','convBC', outIntui,ga_x_mod,csv_Optim, paras_nets{:})

% For AD results, may also use
ShapeDist('AD','convBC', outIntui,dlin,csv_Optim, paras_nets{:})

