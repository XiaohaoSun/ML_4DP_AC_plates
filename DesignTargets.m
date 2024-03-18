
% =========================================================================
% Heading and network
% =========================================================================
addpath(genpath('InverseDesign1'));
addpath('mat_EA_designs');
addpath(genpath('functions'));
addpath('SurfScan');
addpath('netsForward');

% forwardNet
load mat_netD2s_f5convBCxyzEnd_z_12epoDrop.mat
load mat_netD2s_f5convBCxyzEnd_xy_12epoDrop.mat

netxy = netD2s_f5xyz_xy;
netz = netD2s_f5xyz_z;
num_per_voxel = 3;
nx = 15; ny = 15;

global outIntui 
% outIntui is normalized.
% outIntui takes voxel-corner points, 16x16x3.

%% FE-derived target shapes, Figure 5, row 3

flag_flipX = 0;

% --- subR1 ---
% s = rng; 
load('mat_rng_subR1.mat'); rng(s); seed = randi([0,1],5,5,2); temp1 = sub2full(seed,[15,15,2]);
insIntui = [temp1(:,:,1);temp1(:,:,2)];

% FEA generated csv files
csvFile = 'outs-Nocut2_subR1.csv'; frame = 5;
outIntui = readFEA_convBC_frame(csvFile,frame,num_per_voxel); % voxel-corner points
outIntui = (outIntui-cat(3,mux,muy,muz))./cat(3,sigx,sigy,sigz); % normalization

%% FE-derived target shapes, Figure 5, row 2

flag_flipX = 0;

% --- subR2 ---
load('mat_rng_subR2.mat'); rng(s); seed = randi([0,1],3,3,2); temp1 = sub2full(seed,[15,15,2]);
insIntui = [temp1(:,:,1);temp1(:,:,2)];

% FEA generated csv files
csvFile = 'outs-Nocut2_subR2.csv'; frame = 5;
outIntui = readFEA_convBC_frame(csvFile,frame,num_per_voxel); % voxel-corner points
outIntui = (outIntui-cat(3,mux,muy,muz))./cat(3,sigx,sigy,sigz); % normalization

%% FE-derived target shapes, Figure 5, row 1

flag_flipX = 1;

% --- subI4 ---
seed = zeros(5,5,2); seed(1:2,1:2,1)=1; seed(3:5,3:5,1)=1; 
temp1 = sub2full(seed,[15,15,2]);
insIntui = [temp1(:,:,1);temp1(:,:,2)]; % I4

% FEA generated csv files
csvFile = 'outs-Nocut2_subI4.csv'; frame = 5;
outIntui = readFEA_convBC_frame(csvFile,frame,num_per_voxel); % voxel-corner points
if flag_flipX == 1
    % Flip x, 
    insIntuiTemp = cat(3,insIntui(1:ny,:),insIntui(ny+1:end,:));
    [dataX,dataY,dataZ,insIntuiTemp] = data_flipX(outIntui(:,:,1),outIntui(:,:,2),outIntui(:,:,3),insIntuiTemp);
    insIntui = [insIntuiTemp(:,:,1);insIntuiTemp(:,:,2)];
    [dataX,dataY,dataZ] = data_convertBC1n(dataX,dataY,dataZ,2);
    outIntui = cat(3,dataX,dataY,dataZ);
end
outIntui = (outIntui-cat(3,mux,muy,muz))./cat(3,sigx,sigy,sigz); % normalization


%% Visualize

% Visualize design 
visualizeDesign(insIntui);

% Visualize target shape 
colors = tab20;
surf(outIntui(:,:,1)*sigx+mux,outIntui(:,:,2)*sigy+muy,outIntui(:,:,3)*sigz+muz, ...
    'FaceColor',colors(6,:),'FaceAlpha',1,'EdgeColor','none'); 
axis equal; axis off; view(20,45);
% light('Position',[0,-1,0],'Style','infinite');
camlight;
material([0.4 0.5 0.7]); 
lighting gouraud;

