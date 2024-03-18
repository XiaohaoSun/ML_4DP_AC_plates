% EA subdomain design

% Last-step variables needed
% - ga_x                Current optimal solution to start with.
% - outIntui            Target shape.
% - out_Optim_ML_norm   Current optimal shape to identify subdomain.

global outIntui
% outIntui is normalized.
% outIntui takes voxel-corner points, 16x16x3.

%% Skip this section if last-step variables needed are available

% Read results and update ga_x to be the current optimal solution.
open result_check.m
% run 'Read results' section of 'result_check.m'
ga_x = ga_x_mod;

% If ga_x is available, run the following:
out_Optim_ML_norm = forwardNet3(ga_x,netxy,netz,mu);

%% Identify subdomain
% so far, ga_x is the optimal design by global_EA
tol = 25; % Control the subdomain size. 25 used in most cases.
subspace = find_deviation(out_Optim_ML_norm, outIntui, tol);
subspace = conv2(ones(2,2), subspace); subspace = subspace(2:end-1,2:end-1);
subspaceInd = find(cat(3,subspace,subspace));

% Visualize identified subdomain
%{
visualizeSubdomain(subspace)
%}

%% run GA based on forwardNet

ga_popsize = 2000;
maxGen = 5; % 5 or 10

% initial population, partly from ga_x, partly randomized
gaSub_my = ones(400,1)*ga_x(subspaceInd);
gaSub_X0 = [gaSub_my; randi([0,1],ga_popsize-size(gaSub_my,1),length(subspaceInd))]; 

gaSub_nvars = size(gaSub_X0,2); % variable number of an individual
ga_lb = zeros(1,2*nx*ny);
ga_ub = ones(1,2*nx*ny);
IntCon = 1:gaSub_nvars;

% Add custom plot function, population size and population starting point
gaopts = optimoptions('ga','PlotFcn',{"gaplotbestf"},...
                     'PopulationSize',ga_popsize,...
                     'InitialPopulation',gaSub_X0,...
                     'MaxStallGenerations',maxGen,'MaxGenerations',maxGen,...
                     'UseVectorized',true);
% "gaplotstopping" may also be plotted

fprintf('EA started...\n') 
tic
% run Genetic Algorithm
[gaSub_x,gaSub_fval,exitflag,gaSub_out,gaSub_pops] = ...
    ga(@(x)MLfun3sub(x,netD2s_f5xyz_xy,netD2s_f5xyz_z,mu, ga_x, subspaceInd),...
    gaSub_nvars, [],[],[],[], ga_lb,ga_ub,[],IntCon,gaopts);
fprintf('EA finished in %f seconds!\n', toc)

%% save/load EA results

% save results
result_path = fullfile('mat_designs', 'test.mat');
save(result_path,'ga_x','subspace','subspaceInd','gaSub_x','gaSub_pops','gaSub_out');

% load and check results
open result_check.m

%% functions

% ML_shape_prediction for GA
function rmse = MLfun3sub(gaSub_x,netD2Sxy,netD2Sz,mu, ga_x, subspaceInd)
% Support vectorized calculation; mu is for inputs (0 or 1)
global outIntui

ga_temp = ones(size(gaSub_x,1),1)*ga_x;
ga_temp(:,subspaceInd) = gaSub_x;
ga_ins = (reshape(ga_temp',15,15,2,[])-mu)/0.5;

% ML prediction of ins
xyPred = predict(netD2Sxy,ga_ins,'MiniBatchSize',1);
zPred = predict(netD2Sz,ga_ins,'MiniBatchSize',1);
xyzPred = cat(3,xyPred,zPred);

% % back to oriBC?
% [x,y,z,temp] = data_convertBC1_back(xyPred(:,:,1),xyPred(:,:,2),zPred,1);
% xyzPred = cat(3,x,y,z);

rmse = sqrt(mean((xyzPred-outIntui).^2,[1,2,3])); % rmse

% rmse = sqrt(mean(((xyzPred-outIntui).^2)./sqrt(sqrt(1:16)'*sqrt(1:16)),[1,2,3])); % weighted rmse 1
% rmse = sqrt(mean(((xyzPred-outIntui).^2)./sqrt((1:16)'+(1:16)-1),[1,2,3])); % weighted rmse 2
% rmse = sqrt(mean(((xyzPred-outIntui).^2)./((1:16)'+(1:16)-1),[1,2,3])); % weighted rmse 3

rmse = rmse(:); % convert 1x1x1xN to Nx1 vector
end

