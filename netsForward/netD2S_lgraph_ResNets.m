% The code is to build the ResNet models.
% The built ResNet models are already saved, i.e.,
% for z,
% mat_netD2s_lgraph_RN11_singlecoord.mat
% for x,y,
% mat_netD2s_lgraph_RN11_twocoords

% input: 15,15,2
% output 16,16,3

%% Build my ResNet: RN11 (51 conv layers, 24 residual blocks)

% 1. start by modifying ResNet layer
lgraphtemp = resnetLayers([15 15 2],2, BottleneckType="none",...
    InitialStride=1, InitialFilterSize=3, InitialNumFilters=32, InitialPoolingLayer='none', ...
    StackDepth = [13],NumFilters=[64]);
lgraphtemp = removeLayers(lgraphtemp,{'fc','softmax','output'}); % remove some layers for preparation
% analyzeNetwork(lgraphtemp)

% 2. add remaining layers with no ResNet blocks using "replaceLayer"
larray1 = [
    convolution2dLayer(3,128,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,128,'Stride',1,'Padding','same')
    batchNormalizationLayer
    additionLayer(2,'Name','add_block0')
    reluLayer

    convolution2dLayer(3,256,'Stride',1,'Padding',[2 1 2 1])
    batchNormalizationLayer
    reluLayer

    convolution2dLayer(3,128,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,128,'Stride',1,'Padding','same')
    batchNormalizationLayer
    additionLayer(2,'Name','add_block1')
    reluLayer

    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    additionLayer(2,'Name','add_block2')
    reluLayer

    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    additionLayer(2,'Name','add_block3')
    reluLayer

    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    additionLayer(2,'Name','add_block4')
    reluLayer

    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    additionLayer(2,'Name','add_block5')
    reluLayer

    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    additionLayer(2,'Name','add_block6')
    reluLayer

    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    additionLayer(2,'Name','add_block7')
    reluLayer

    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    additionLayer(2,'Name','add_block8')
    reluLayer

    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    additionLayer(2,'Name','add_block9')
    reluLayer

    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    additionLayer(2,'Name','add_block10')
    reluLayer

    convolution2dLayer(3,2,'Stride',1,'Padding','same')
    regressionLayer];

lgraph = replaceLayer(lgraphtemp,'gap',larray1);

% 3. make ResNet blocks with added layers
% 3.0 make my resnet block 0 before 15x15->16x16 layer
nameLastRelu = lgraphtemp.Layers(end-1).Name;
skipConv = [
    convolution2dLayer(1,128,'Stride',1,'Padding','same','Name','skipConv0')
    batchNormalizationLayer('Name','skipBN0')];
lgraph = addLayers(lgraph,skipConv);
lgraph = connectLayers(lgraph,nameLastRelu,'skipConv0');
lgraph = connectLayers(lgraph,'skipBN0','add_block0/in2');

% 3.1 make my resnet block 1
skipConv = [
    convolution2dLayer(1,128,'Stride',1,'Padding','same','Name','skipConv1')
    batchNormalizationLayer('Name','skipBN1')];
lgraph = addLayers(lgraph,skipConv);
lgraph = connectLayers(lgraph,'relu_3','skipConv1');
lgraph = connectLayers(lgraph,'skipBN1','add_block1/in2');

% 3.2 make my resnet block 2
skipConv = [
    convolution2dLayer(1,64,'Stride',1,'Padding','same','Name','skipConv2')
    batchNormalizationLayer('Name','skipBN2')];
lgraph = addLayers(lgraph,skipConv);
lgraph = connectLayers(lgraph,'relu_5','skipConv2');
lgraph = connectLayers(lgraph,'skipBN2','add_block2/in2');

% 3.3 make my resnet block 3
lgraph = connectLayers(lgraph,'relu_7','add_block3/in2');

% 3.4 make my resnet block 4
lgraph = connectLayers(lgraph,'relu_9','add_block4/in2');

% make my resnet block 5
lgraph = connectLayers(lgraph,'relu_11','add_block5/in2');

% make my resnet block 6
lgraph = connectLayers(lgraph,'relu_13','add_block6/in2');

% make my resnet block 7
lgraph = connectLayers(lgraph,'relu_15','add_block7/in2');

% make my resnet block 8
lgraph = connectLayers(lgraph,'relu_17','add_block8/in2');

% make my resnet block 9
lgraph = connectLayers(lgraph,'relu_19','add_block9/in2');

% make my resnet block 10
lgraph = connectLayers(lgraph,'relu_21','add_block10/in2');

analyzeNetwork(lgraph)

