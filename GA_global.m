% EA global design

% Variables needed
% - outIntui            Target shape.

global outIntui 
% outIntui is normalized.
% outIntui takes voxel-corner points, 16x16x3.

%{
open DesignTargets.m
%}

%% run GA based on forwardNet

ga_popsize = 2000;
maxGen = 100;
maxStallGen = 50;

% Two types of initial populations
% 1. initial population fully random
ga_X0 = randi([0,1],ga_popsize,2*nx*ny); 

% % 2. initial population by inverseNet
% ga_X0 = [ga_my;randi([0,1],ga_popsize-size(ga_my,1),2*nx*ny)]; 

ga_nvars = size(ga_X0,2); % variable number of a individual
ga_lb = zeros(1,2*nx*ny);
ga_ub = ones(1,2*nx*ny);
IntCon = 1:ga_nvars;

% Add custom plot function, population size and population starting point
gaopts = optimoptions('ga','PlotFcn',{"gaplotbestf"},...
                     'PopulationSize',ga_popsize,...
                     'InitialPopulation',ga_X0,...
                     'MaxStallGenerations',maxStallGen,'MaxGenerations',maxGen,...
                     'UseVectorized',false);
% "gaplotstopping" may also be plotted

fprintf('EA started...\n') 
tic
% run Genetic Algorithm
[ga_x,ga_fval,exitflag,ga_out,ga_pops] = ...
    ga(@(x)MLfun3(x,netD2s_f5xyz_xy,netD2s_f5xyz_z,mu), ga_nvars, [],[],[],[],...
    ga_lb,ga_ub,[],IntCon,gaopts);
fprintf('EA finished in %f seconds!\n', toc)

%% save/load EA results

% save results
result_path = fullfile('mat_designs', 'test.mat');
save(result_path,'ga_x','subspace','subspaceInd','gaSub_x','gaSub_pops','gaSub_out');

% load and check results
open result_check.m

%% functions

% ML_shape_prediction for GA
function rmse = MLfun2(ga_x,netD2S,netD2Sz,mu)
% Support vectorized calculation; mu is for inputs (0 or 1)

global outIntui
ga_ins = reshape(ga_x',15,15,2,[]);
ga_ins = (ga_ins-mu)/0.5; % normalization for ins (mu is not exactly 0.5)

% ML prediction of ins
xyzPred = predict(netD2S,ga_ins,'MiniBatchSize',250);
zPred = predict(netD2Sz,ga_ins,'MiniBatchSize',250);

xyzPred = cat(3,xyzPred(:,:,1:2),zPred);
rmse = sqrt(mean((xyzPred-outIntui).^2,[1,2,3])); % rmse
rmse = rmse(:); % convert 1x1x1xN to Nx1 vector
end


% ML_shape_prediction for GA
function rmse = MLfun3(ga_x,netD2Sxy,netD2Sz,mu)
% Support vectorized calculation; mu is for inputs (0 or 1)

global outIntui
ga_ins = (reshape(ga_x',15,15,2,[])-mu)/0.5;

% ML prediction of ins
xyPred = predict(netD2Sxy,ga_ins,'MiniBatchSize',250);
zPred = predict(netD2Sz,ga_ins,'MiniBatchSize',250);
xyzPred = cat(3,xyPred,zPred);

% % back to oriBC?
% [x,y,z,temp] = data_convertBC1_back(xyPred(:,:,1),xyPred(:,:,2),zPred,1);
% xyzPred = cat(3,x,y,z);


% rmse = sqrt(mean((xyzPred-outIntui).^2,[1,2,3])); % rmse

% rmse = sqrt(mean(((xyzPred-outIntui).^2)./sqrt(sqrt(1:16)'*sqrt(1:16)),[1,2,3])); % weighted rmse 1

rmse = sqrt(mean(((xyzPred-outIntui).^2)./sqrt((1:16)'+(1:16)-1),[1,2,3])); % weighted rmse 2

% rmse = sqrt(mean(((xyzPred-outIntui).^2)./((1:16)'+(1:16)-1),[1,2,3])); % weighted rmse 3


rmse = rmse(:); % convert 1x1x1xN to Nx1 vector
end


