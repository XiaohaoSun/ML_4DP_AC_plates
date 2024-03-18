% GD design, gradient by AD.

% Variables needed
% - outIntui            Target shape.

global outIntui 
% outIntui is normalized.
% outIntui takes voxel-corner points, 16x16x3.

%{
open DesignTargets.m
%}

%% Optimization, custom trainning loop
% =========================================================================
% Initial patterns, four types
% =========================================================================

% all passive initilization
in0 = (reshape(zeros(1,2*nx*ny),nx,nx,2)-mu)/0.5;

% all active initilization
in0 = (reshape(ones(1,2*nx*ny),nx,nx,2)-mu)/0.5;

% random initilization
rng('default'); % put here to make results reproducible; may be deleted.
in0 = (randi([0,1],[nx,nx,2])-mu)/0.5;

% neutral initialization
in0 = zeros(nx,nx,2);

% =========================================================================


% Convert to dlnetwork objects
dlin = dlarray(in0,'SSCB'); 
netxy = netD2s_f5xyz_xy;
netz = netD2s_f5xyz_z;
[dlnetxy, state0_netxy] = convert_netDAG_to_dlnet(netxy);
[dlnetz, state0_netz] = convert_netDAG_to_dlnet(netz);
dloutTarget = dlarray(outIntui,'SSCB');

%% Optimization, custom trainning loop

executionEnvironment = "auto"; 

if canUseGPU
    dlin = gpuArray(dlin);
    dloutTarget = gpuArray(dloutTarget);
end

clearCache(accfun); % Added on Apr.13, previous results have not used it.

% accfun = dlaccelerate(@objectiveAndGradient1);
% accfun = dlaccelerate(@objectiveAndGradient2);
accfun = dlaccelerate(@objectiveAndGradient3);


plots = "training-progress";

% Initilize the training progress plot
if plots == "training-progress"
    figure
    lineLossTrain = animatedline('Color',[0.85 0.325 0.098]);
    ylim([0 inf])
    xlabel("Iteration")
    ylabel("Loss")
    grid on
end

numEpochs = 100; 
miniBatchSize = 1;

% Initialize the average gradients and squared average gradients for the ADAM solver
initialLearnRate = 10*(1e-3);
learnRateDropFactor = 1; % sqrt(0.2);

% Initialize some parameters for the ADAM solver
gradDecay = 0.9;
sqGradDecay = 0.99; % default is 0.999
averageGrad = [];
averageSqGrad = [];

numIterationsPerEpoch = 10;

iteration = 0;
start = tic;
loss_history = [];

% Loop over epochs.
for epoch = 1:numEpochs
    
    % Determine learning rate.
    if epoch == 1
        learnRate = initialLearnRate;
    elseif rem(epoch-1,50) == 0
        learnRate = learnRate*learnRateDropFactor;
    end
    
    % Loop over mini-batches.
    for i = 1:numIterationsPerEpoch
        iteration = iteration + 1;
        
        
        % Evaluate the model gradients, state, and loss using dlfeval and the
        % modelGradients function and update the network state.
%         [loss,gradients] = dlfeval(@objectiveAndGradient3,dlin,dlnetxy,dlnetz,dloutTarget);
        [loss,gradients] = dlfeval(accfun,dlin,dlnetxy,dlnetz,dloutTarget);
        
        % Update the network parameters using the SGDM optimizer.
        [dlin, averageGrad,averageSqGrad] = adamupdate(dlin, gradients, ...
            averageGrad, averageSqGrad, iteration, learnRate, gradDecay, sqGradDecay);
        
        % Display the training progress.
        if plots == "training-progress"
            loss = double(gather(extractdata(loss)));
            addpoints(lineLossTrain,iteration,loss)
            D = duration(0,0,toc(start),'Format','hh:mm:ss');
            title("Epoch: " + epoch + ", Elapsed: " + string(D) + ", Loss: " + loss)
            drawnow
        end
        % store loss history
        loss_history(epoch) = loss;
    end
end

% % Check
% accfun

%% save/load AD results

% save results
result_path = fullfile('mat_designs', 'test.mat');
save(result_path,'dlin','dlnetxy','dlnetz','outIntui');

% load and check results
open result_check.m

%% functions
function [loss, grad] = objectiveAndGradient1(dlin,dlnetxy,dlnetz,dloutTarget)
    dlxyPred = predict(dlnetxy,dlin);
    dlzPred = predict(dlnetz,dlin);
    loss1 = sum((dloutTarget - cat(3,dlxyPred,dlzPred)).^2,'all')/768;
%     loss1 = sum( ((dloutTarget - [dlxPred;dlyPred]).^2)./(ones(2,1)*sqrt(1:24)), 'all')...
%     /sum(ones(2,1)*sqrt(1:24),'all'); % distance-weighted?
    loss2 = mean(((2-2*0.4993336-dlin).*(2*0.4993336+dlin)).^2,'all');
    loss = 0.9*loss1 + 0.2*loss2;
    grad = dlgradient(loss,dlin);
end

function [loss, grad] = objectiveAndGradient2(dlin0,dlnetxy,dlnetz,dloutTarget)
    dlin = tanh(dlin0);
    dlxyPred = predict(dlnetxy,dlin);
    dlzPred = predict(dlnetz,dlin);
    loss1 = sum((dloutTarget - cat(3,dlxyPred,dlzPred)).^2,'all')/768;
%     loss1 = sum( ((dloutTarget - [dlxPred;dlyPred]).^2)./(ones(2,1)*sqrt(1:24)), 'all')...
%     /sum(ones(2,1)*sqrt(1:24),'all'); % distance-weighted?
    loss2 = mean((1-dlin.^2),'all');
    loss = 0.9*loss1 + 0.1*loss2;
    grad = dlgradient(loss,dlin);
end

function [loss, grad] = objectiveAndGradient3(dlin0,dlnetxy,dlnetz,dloutTarget)
    dlin = tanh(dlin0);
    dlin = dlin + 1 - 2*0.4993336; % (-1,1) -> (-2*mu, 2-2*mu)
    dlxyPred = predict(dlnetxy,dlin);
    dlzPred = predict(dlnetz,dlin);
    loss1 = sum((dloutTarget - cat(3,dlxyPred,dlzPred)).^2,'all')/768;
%     loss1 = sum( ((dloutTarget - [dlxPred;dlyPred]).^2)./(ones(2,1)*sqrt(1:24)), 'all')...
%     /sum(ones(2,1)*sqrt(1:24),'all'); % distance-weighted?
    loss2 = mean((2-2*0.4993336-dlin).*(2*0.4993336+dlin),'all');
    loss = 0.9*loss1 + 0.1*loss2;
    grad = dlgradient(loss,dlin);
end


function [dlnet, state_dlnet] = convert_netDAG_to_dlnet(net)
lgraph_net = layerGraph(net); 
lgraph_net = removeLayers(lgraph_net,'regressionoutput');
dlnet = dlnetwork(lgraph_net);
state_dlnet = dlnet.State;
end


function dlinsPred = modelPredictions1(dlnet,YTest)
dlYTest = dlarray(reshape(cell2mat(YTest'),size(YTest{1},1),size(YTest{1},2),numel(YTest)),'CTB');
dlinsPred = predict(dlnet,dlYTest);
end

function dlinsPred = modelPredictions(dlnet,YTest)
dlYTest = dlarray(reshape(cell2mat(YTest'),size(YTest{1},1),size(YTest{1},2),numel(YTest)),'CTB');
dlinsPred = predict(dlnet,dlYTest);
% extractdata(dlinsPred)
end

% ML forward shape prediction
function ga_xyzPred = forwardNet3(ga_x,netD2Sxy,netD2Sz,mu)
ga_ins = reshape(ga_x,[15,15,2]);
ga_ins = (ga_ins-mu)/0.5; % normalization for ins (mu is not exactly 0.5)
xyPred = predict(netD2Sxy,ga_ins,'MiniBatchSize',1);
zPred = predict(netD2Sz,ga_ins,'MiniBatchSize',1);
ga_xyzPred = cat(3,xyPred,zPred);
end
