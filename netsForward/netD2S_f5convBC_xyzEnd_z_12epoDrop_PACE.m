% design to shape, CNN
% adam
clear all
close all
clc

addpath(genpath('functions'));
addpath(genpath('dataset'));

%% dataset
% frame 5 only

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

% Shuffle the data, must be excuted following the data loading
% seqNew = randperm(size(inBinary,4),size(inBinary,4));
load('data_shuffle_seq_220704_900k.mat');
outx = outx(:,:,:,seqNew);
outy = outy(:,:,:,seqNew);
outz = outz(:,:,:,seqNew);
inBinary = inBinary(:,:,:,seqNew);
% % ----------------------------------------

% Take voxel's feature coordinates
nx = 15; ny = 15; 
num_per_voxel_x = (size(outx,1)-1)/nx;
num_per_voxel_y = (size(outx,2)-1)/ny;

% Corner-point
outz = outz((1):num_per_voxel_x:end,(1):num_per_voxel_y:end,:,:);
outy = outy((1):num_per_voxel_x:end,(1):num_per_voxel_y:end,:,:);
outx = outx((1):num_per_voxel_x:end,(1):num_per_voxel_y:end,:,:);

% standardize data for input
mu = mean(inBinary(:));
sig = 0.5; % std(inBinary(:),0,2);
inBinary = (inBinary - mu)/sig;

% standardize data for output x y z
[mux,sigx,outx] = NormData(outx);
[muy,sigy,outy] = NormData(outy);
[muz,sigz,outz] = NormData(outz);

% We choose: 90% training, 10% validation data
numData = size(outx,4);
miniBatchSize  = 512; 
numTrainSet = floor(round(numData*0.9)/miniBatchSize)*miniBatchSize;

clear outx outy

% insTrain = inBinary(:,:,:,1:numTrainSet); % inputs Train
% f5_XYZTrain = outz(:,:,:,1:numTrainSet); %
insValid = inBinary(:,:,:,numTrainSet+1:end); % inputs Valid or Test
f5_XYZValidz = outz(:,:,:,numTrainSet+1:end); % 

fprintf('Data preparation completed. This step takes %f seconds!\n', toc) 

%% Build network
% layers constructed in local machine since PACE does not support R2022a.
% On PACE, we load the built lgraph
load('mat_netD2s_lgraph_RN11_singlecoord.mat');

validationFrequency = floor(numTrainSet/miniBatchSize);
% sgdm or adam
option1 = trainingOptions('adam', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',80, ...
    'InitialLearnRate',1e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',sqrt(0.2), ...
    'LearnRateDropPeriod',12, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{insValid,f5_XYZValidz}, ...
    'ValidationFrequency',validationFrequency, ...
    'Plots','none', ... or training-progress,
    'Verbose',true, ... VerboseFrequency 50 by default
    'CheckpointPath','Checkpoints'); % one per epoch by default

% Plots none, Verbose true for non-interactive mode.
% Plots training-progress, Verbose false for GUI mode.

%% Training
tic
[netD2s_f5xyz_z,info] = trainNetwork(inBinary(:,:,:,1:numTrainSet),...
    outz(:,:,:,1:numTrainSet),lgraph,option1);
fprintf('Net training completed. This step takes %f seconds!\n', toc) 

tic
save('mat_netD2s_f5convBCxyzEnd_RN22z_51conv_12epoDrop.mat','netD2s_f5xyz_z','info',...
    'mu','mux','muy','muz','sig','sigx','sigy','sigz',...
    'numData','numTrainSet','miniBatchSize',...
    'insValid','f5_XYZValidz','-v7.3');
fprintf('Resutls saved within %f seconds!\n', toc)

