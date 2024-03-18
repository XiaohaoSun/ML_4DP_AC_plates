% design to shape, GCN
% All on convBC1
% All on xyzEnd (corner-point)
% All on z only

% Here, SC4 has 4 SCs, similar to GCN-Res in Kipf paper or DeepGCN paper.
% Each SC skips multiple (4 or 5 in this case) layers, 
% and Z_n is outside ReLU in the form Z_m = relu(...) + Z_n.

clear all
close all
clc

addpath(genpath('functions'));
addpath(genpath('dataset'));

%% dataset

tic

load('Dataset_frame5_random_design_single.mat'); 
load('Dataset_frame5_island_design_single.mat'); % not finish extracting

outx = cat(4,outx220704rand,outx220704);
outy = cat(4,outy220704rand,outy220704);
outz = cat(4,outz220704rand,outz220704);
inBinary = cat(4,inBinary220704rand,inBinary220704);
clear inBinary220704 inBinary220704rand outx220704 outx220704rand
clear outy220704 outy220704rand outz220704 outz220704rand

[outx,outy,outz,inBinary] = data_augmentation(outx,outy,outz,inBinary);
[outx,outy,outz,inBinary] = data_convertBC1(outx,outy,outz,inBinary);

% Shuffle the data
% seqNew = randperm(size(inBinary,4),size(inBinary,4));
load('data_shuffle_seq_220704_900k.mat');
outx = outx(:,:,:,seqNew);
outy = outy(:,:,:,seqNew);
outz = outz(:,:,:,seqNew);
inBinary = inBinary(:,:,:,seqNew);

% Take voxel's feature coordinates
nx = 15; ny = 15;
num_per_voxel_x = (size(outz,1)-1)/nx;
num_per_voxel_y = (size(outz,2)-1)/ny;

% Corner-point
outx = outx((1):num_per_voxel_x:end,(1):num_per_voxel_y:end,:,:);
outy = outy((1):num_per_voxel_x:end,(1):num_per_voxel_y:end,:,:);
outz = outz((1):num_per_voxel_x:end,(1):num_per_voxel_y:end,:,:);

% standardize data for input
mu = mean(inBinary(:)); % not directly use 0.5
sig = 0.5; % std(inBinary(:),0,2);
inBinary = (inBinary - mu)/sig;

% standardize data for output x y z
[mux,sigx,outx] = NormData(outx);
[muy,sigy,outy] = NormData(outy);
[muz,sigz,outz] = NormData(outz);

% (Be careful about standardizing data)

% define a constant adjacency matrix, 15x15 nodes
tempA= reshape(1:(ny*nx), ny, nx);
G = graph();
for ix = 1:(nx-1)
    G = addedge(G, tempA(:,ix), tempA(:,ix+1));
end 
for iy = 1:(ny-1)
    G = addedge(G, tempA(iy,:), tempA(iy+1,:));
end
adjMatrix = adjacency(G);

% Initial node locations
[outx0, outy0] = meshgrid(1:15); % not actual x0,y0, but indices.
outx0y0 = cat(3, outx0, outy0);
[mux0y0, sigx0y0, outx0y0] = NormData(outx0y0);

% use outz later, and clear outx and outy
clear outx outy

% We choose: 90% training, 10% validation data
numData = size(outz,4); % 899744 for random + island
miniBatchSize  = 512; 
numTrainSet = floor(round(numData*0.9)/miniBatchSize)*miniBatchSize;

% Make array datastores, for the use of mini-batches of data.
dsTrain_ins = arrayDatastore(inBinary(:,:,:,1:numTrainSet), "IterationDimension",4);
dsTrain_outs = arrayDatastore(outz(:,:,:,1:numTrainSet), "IterationDimension",4);
dsTrain = combine(dsTrain_ins, dsTrain_outs);

% For huge validation set, do not process whole data directly,
% ---but process inputs in mini-batches
insValid = inBinary(:,:,:,numTrainSet+1:numData);
dsValid = arrayDatastore(insValid, "IterationDimension",4); % only ins in datastore
dsValid.ReadSize = miniBatchSize;
% ---and process outputs in whole.
fprintf('size of outz is %d, %d, %d, %d.\n', size(outz));
g_f5_XYZValidz = myPreprocessOuts(outz(:,:,:,numTrainSet+1:numData));

fprintf('Data preparation completed. This step takes %f seconds!\n', toc) 

%% define GCN and prepare for training
% (later, may explore other GCN architectures, or GATs)
parameters = struct;
numHiddenFeatureMaps = 64;

numInputFeatures = size(inBinary, 3) + size(outx0y0, 3);
numOuputFeagures = size(outz, 3);

parameters.mult1.Weights = myInitializeGlorot( ...
    [numInputFeatures, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numInputFeatures, 'double');

parameters.mult2.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult3.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult4.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult5.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult6.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult7.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult8.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult9.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult10.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult11.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult12.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult13.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult14.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult15.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult16.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult17.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult18.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult19.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult20.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numHiddenFeatureMaps], ...
    numHiddenFeatureMaps, numHiddenFeatureMaps, 'double');

parameters.mult21.Weights = myInitializeGlorot( ...
    [numHiddenFeatureMaps, numOuputFeagures], ...
    numOuputFeagures, numHiddenFeatureMaps, 'double');


% TRAINING OPTIONS
numEpochs = 20;
initialLearnRate = 0.01;
LRDropFactor = sqrt(0.2);
LRDropPeriod = 3;
validationFrequency = floor(numTrainSet/miniBatchSize);
verboseFrequency = 50; % change to 50 for real training on pace

% PREPARE TRAINING MINIBATCHES
mbq = minibatchqueue(dsTrain, 3, ...
    MiniBatchSize=miniBatchSize, ...
    PartialMiniBatch="discard", ...
    MiniBatchFcn=@(ins, outs) myPreprocessMiniBatch(adjMatrix, outx0y0, ins, outs), ...
    OutputCast="double", ...
    OutputAsDlarray=[0 1 0], ...
    OutputEnvironment = ["cpu", "auto", "cpu"]);


%% TRAIN MODEL

% INITIALIZE SOME ADAM PARAMETERS
trailingAvg = [];
trailingAvgSq = [];

% % dlaccelerate does not speed the training up in this case.
% accfun = dlaccelerate(@modelLoss);
% clearCache(accfun);
% accfun = dlaccelerate(@modelLoss);

% INITIALIZE TRAINING INFORMATION 
info_GCN = struct;
numIterations = floor(numTrainSet/miniBatchSize)*numEpochs;
info_GCN.TrainingLoss = NaN(numIterations, 1);
info_GCN.ValidationLoss = info_GCN.TrainingLoss;
info_GCN.BaseLearnRate = info_GCN.TrainingLoss;
info_GCN.totalTime = NaN;

% START TRAINING
fprintf('Start training the ML model.\n')
fprintf('|=============================================================================|\n');
fprintf('|  Epoch  | Iteration | Time Elapsed | Mini-batch | Validation |  Base Learn  |\n')
fprintf('|         |           |  (hh:mm:ss)  |    Loss    |    Loss    |     Rate     |\n')
fprintf('|=============================================================================|\n');

iteration = 0;
start = tic;
for epoch = 1:numEpochs

    % LearnRate drop
    if epoch == 1
        learnRate = initialLearnRate;
    elseif rem(epoch-1, LRDropPeriod) == 0
        learnRate = learnRate * LRDropFactor;
    end

    % Shuffle data
    shuffle(mbq);

    while hasdata(mbq)
        iteration = iteration + 1;

        % Read mini-batches of data
        [g_adjTrain, g_insTrain, g_f5_XYZTrain] = next(mbq);

        % Evaluate the model loss and gradients.
        [loss,gradients] = dlfeval(@modelLoss,parameters,g_insTrain,g_adjTrain,g_f5_XYZTrain);
%         [loss,gradients] = dlfeval(accfun,parameters,insTrain,adjTrain,f5_XYZTrain);
    
        % Update the network parameters using the Adam optimizer.
        [parameters,trailingAvg,trailingAvgSq] = adamupdate(parameters,gradients, ...
            trailingAvg,trailingAvgSq,epoch,learnRate);

        % Record infos
        info_GCN.TrainingLoss(iteration) = loss;
        info_GCN.BaseLearnRate(iteration) = learnRate;
    
        % Display the validation metrics.
        if iteration == 1 || mod(iteration,validationFrequency) == 0
            % make validation predictions on mini-batches, using "read".
            g_f5_XYZPredz = [];
            reset(dsValid);
            while hasdata(dsValid)
                data = read(dsValid);
                % the read data are cell arrays N x M, 
                % where N is the number of observations in a mini-batch,
                % and M (here, = 1) is the number of datastores in dsValid. 
                [mb_g_adjValid, mb_g_insValid] = ...
                    myPreprocessMiniBatch(adjMatrix, outx0y0, data(:,1));
                % predict the mini-batch and store
                g_f5_XYZPredz = [g_f5_XYZPredz; model(parameters,mb_g_insValid,mb_g_adjValid)];
            end
            lossValidation = mse(g_f5_XYZPredz,g_f5_XYZValidz,DataFormat="BC");
            info_GCN.ValidationLoss(iteration) = lossValidation;

            % Print validation information on screen
            D = duration(0,0,toc(start),Format="hh:mm:ss");
            fprintf('| %7d | %9d | %12s | %10.4f | %10.4f | %12.4e |\n',...
                epoch, iteration, string(D), loss, lossValidation, learnRate);

        end

        % Print verbose information on screen
        if mod(iteration, verboseFrequency) == 0
            D = duration(0,0,toc(start),Format="hh:mm:ss");
            fprintf('| %7d | %9d | %12s | %10.4f |            | %12.4e |\n',...
                epoch, iteration, string(D), loss, learnRate);
        end

    end

end
fprintf('|=============================================================================|\n');
info_GCN.totalTime = toc(start);
fprintf('Net training completed. This step takes %f seconds!\n', info_GCN.totalTime); 

%% Use PACE to make predictions on validation dataset (which is huge)

tic
g_f5_XYZPredz = [];
reset(dsValid);
while hasdata(dsValid)
    data = read(dsValid);
    % the read data are cell arrays N x M, 
    % where N is the number of observations in a mini-batch,
    % and M (here, = 1) is the number of datastores in dsValid. 
    [mb_g_adjValid, mb_g_insValid] = ...
        myPreprocessMiniBatch(adjMatrix, outx0y0, data(:,1));
    % predict the mini-batch and store
    g_f5_XYZPredz = [g_f5_XYZPredz; model(parameters,mb_g_insValid,mb_g_adjValid)];
end
fprintf('Final prediction of validation dataset completed, obtaining g_f5_XYZPredz. This step takes %f seconds!\n', toc); 

%% Save model and some parameters

tic
save('mat_GCN3z_Acht21x64_SC4.mat','parameters','info_GCN',...
    'mu','muz','sig','sigz','mux0y0','sigx0y0','outx0y0','adjMatrix',...
    'numData','numTrainSet','miniBatchSize',...
    'insValid', 'g_f5_XYZPredz','g_f5_XYZValidz','-v7.3');
fprintf('Resutls saved within %f seconds!\n', toc)

% On Pace training, make sure
% ---always use "double" in parameters initialization and mbq call.
% ---miniBatchSize = 512
% ---verboseFrequency = 50
% ---plotFlag = false
% ---Store insValid, g_f5_XYZPredz, g_f5_XYZValidz for later use.
%    Then, no need to store g_adjValid and g_insValid.

%% FUNCTIONS

function [g_adjacency, g_features, g_outs] = myPreprocessMiniBatch(adjMatrix, x0y0, inBinary, outs)
    % Process mini-batches of data, used with datastore or minibatchqueue.
    % Arguments:
    % ---"adjMatrix" and "x0y0" are the same for all N observations.
    %     adjMatrix, 2D sparse adjacency matrix.
    %     x0y0, initial index coordinates.
    % ---"inBinary" and "outs" are N-by-1 cell arrays,
    %     where N is the number of observations in the mini-batch.
    
    % === Process inputs to give "g_adjacency" and "g_features" ===
    % Concatenate data in cell arrays along the 4-th dimision, 
    % obtaining "inBinary" in 4D, 15x15x2xN (not padded to 16x16x2xN).
    inBinary = cat(4, inBinary{:});
    % Then process.
    g_adjacency = sparse([]);
    numData = size(inBinary,4);
    g_features = -1*ones(15*15*numData,size(inBinary,3)+size(x0y0,3)); % features contain inBinary and x0y0. 
    for idx = 1:size(inBinary,4)
        g_adjacency = blkdiag(g_adjacency, adjMatrix);
        g_features( (15*15*(idx-1)+1):(15*15*idx), :) = cat(2, ...
            reshape(inBinary(:,:,:,idx),15*15,[]), reshape(x0y0, 15*15, []) );
    end

    % === Process "outs" data to give "g_outs" ===
    if nargin > 3
        % Concatenate data in cell arrays along the 4-th dimision,
        % obtaining "outs" in 4D, 16x16x(1,2 or 3)xN.
        outs = cat(4, outs{:});
        % Then process.
        outs = outs(2:end,2:end,:,:);
        g_outs = -1*ones(15*15*numData,size(outs,3));
        for idx = 1:size(outs,4)
            g_outs( (15*15*(idx-1)+1):(15*15*idx), : ) = reshape( outs(:,:,:,idx), 15*15, []);
        end
    end

end

function g_outs = myPreprocessOuts(outs)
    % See details in "myPreprocessWhole".
    % select "outs", from 16x16x(1,2 or 3)xN to 15x15x(1,2 or 3)xN
    outs = outs(2:end,2:end,:,:);
    numData = size(outs,4);
    g_outs = -1*ones(15*15*numData,size(outs,3));
    for idx = 1:numData
        g_outs( (15*15*(idx-1)+1):(15*15*idx), : ) = reshape( outs(:,:,:,idx), 15*15, []);
    end
end

function [g_adjacency, g_features] = myPreprocessIns(adjMatrix, x0y0, inBinary)
    % See details in "myPreprocessWhole".
    g_adjacency = sparse([]);
    numData = size(inBinary,4);
    g_features = -1*ones(15*15*numData,size(inBinary,3)+size(x0y0,3)); % features contain inBinary and x0y0.
    for idx = 1:numData
        g_adjacency = blkdiag(g_adjacency, adjMatrix);
        g_features( (15*15*(idx-1)+1):(15*15*idx), :) = cat(2, ...
            reshape(inBinary(:,:,:,idx),15*15,[]), reshape(x0y0, 15*15, []) );
    end
end

function [g_adjacency, g_features, g_outs] = myPreprocessWhole(adjMatrix, x0y0, inBinary, outs)
    % Process whole data: whole inputs and outputs if nargin==4; 
    %                     whole inputs only if nargin==3.
    % Arguments:
    % ---adjMatrix, 2D sparse matrix, same for all N observations.
    % ---x0y0, initial index coordinates, same for all N observations.
    % ---inBinary, 4D data, 15x15x2xN (not padded to 16x16x2xN).
    % ---outs, 4D data, 16x16x(1,2 or 3)xN.
    
    % Process inputs to give "g_adjacency" and "g_features"
    [g_adjacency, g_features] = myPreprocessIns(adjMatrix, x0y0, inBinary);

    % Process outputs ("outs") data to give "g_outs".
    if nargin > 3
        g_outs = myPreprocessOuts(outs);
    end

end


% MODEL FUNCTION
function [Y, Z2, Z1] = model(parameters,X,A)
% Note:
% It seems ReLU is better than ELU.

ANorm = normalizeAdjacency(A);

Z1 = X;

Z2 = ANorm * Z1 * parameters.mult1.Weights;
% Z2 = relu(Z2) + Z1; % This looks similar to skip connection, 
%                     % but Z1's channel size=2, not compatible with Z2.
% Z2 = relu(Z2) + Z1(:,1) - Z1(:,2);
Z2 = relu(Z2);

Z3 = ANorm * Z2 * parameters.mult2.Weights;
Z3 = relu(Z3);

Z4 = ANorm * Z3 * parameters.mult3.Weights;
Z4 = relu(Z4);

Z5 = ANorm * Z4 * parameters.mult4.Weights;
Z5 = relu(Z5);

Z6 = ANorm * Z5 * parameters.mult5.Weights;
Z6 = relu(Z6) + Z2;

Z7 = ANorm * Z6 * parameters.mult6.Weights;
Z7 = relu(Z7);

Z8 = ANorm * Z7 * parameters.mult7.Weights;
Z8 = relu(Z8);

Z9 = ANorm * Z8 * parameters.mult8.Weights;
Z9 = relu(Z9);

Z10 = ANorm * Z9 * parameters.mult9.Weights;
Z10 = relu(Z10);

Z11 = ANorm * Z10 * parameters.mult10.Weights;
Z11 = relu(Z11) + Z6;

Z12 = ANorm * Z11 * parameters.mult11.Weights;
Z12 = relu(Z12);

Z13 = ANorm * Z12 * parameters.mult12.Weights;
Z13 = relu(Z13);

Z14 = ANorm * Z13 * parameters.mult13.Weights;
Z14 = relu(Z14);

Z15 = ANorm * Z14 * parameters.mult14.Weights;
Z15 = relu(Z15);

Z16 = ANorm * Z15 * parameters.mult15.Weights;
Z16 = relu(Z16) + Z11;

Z17 = ANorm * Z16 * parameters.mult16.Weights;
Z17 = relu(Z17);

Z18 = ANorm * Z17 * parameters.mult17.Weights;
Z18 = relu(Z18);

Z19 = ANorm * Z18 * parameters.mult18.Weights;
Z19 = relu(Z19);

Z20 = ANorm * Z19 * parameters.mult19.Weights;
Z20 = relu(Z20);

Z21 = ANorm * Z20 * parameters.mult20.Weights;
Z21 = relu(Z21) + Z16;

Z22 = ANorm * Z21 * parameters.mult21.Weights;

Y = Z22;

end

% MODEL LOSS FUNCTION
function [loss,gradients] = modelLoss(parameters,X,A,T)

Y = model(parameters,X,A);
loss = mse(Y,T,DataFormat="BC");
gradients = dlgradient(loss, parameters);

end

% REDEFINE myInitializeGlorot and normalizeAdjacency here to avoid problems.
function weights = myInitializeGlorot(sz,numOut,numIn,className)

arguments
    sz
    numOut
    numIn
    className = 'single'
end

Z = 2*rand(sz,className) - 1;
bound = sqrt(6 / (numIn + numOut));

weights = bound * Z;
weights = dlarray(weights);

end

function ANorm = normalizeAdjacency(A)

% Add self connections to adjacency matrix.
A = A + speye(size(A));

% Compute inverse square root of degree.
degree = sum(A, 2);
degreeInvSqrt = sparse(sqrt(1./degree));

% Normalize adjacency matrix.
ANorm = diag(degreeInvSqrt) * A * diag(degreeInvSqrt);

end