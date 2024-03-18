function [out_Optim_ML,out_Optim_FE] = ShapeDistP2P(optimType,BCType, ...
    outIntui, in_Optim, csv_Optim, ...
    netxy,netz, mu,mux,muy,muz,sig,sigx,sigy,sigz)
% optimType = 'AD' or 'GA'
% BCType = 'convBC' or 'oriBC'

% target
outIntuiOri = outIntui.*cat(3,sigx,sigy,sigz)+cat(3,mux,muy,muz);

% optim_ML
if strcmp(optimType,'AD')
%     out_Optim_ML_norm = forwardNet3(round_of_dlin(in_Optim)*sig+mu,netxy,netz,mu); % type 2
    out_Optim_ML_norm = forwardNet3(round_of_dlin(in_Optim)*0.5+0.5,netxy,netz,mu); % type 3
elseif strcmp(optimType,'GA')
    % out_Optim_ML_norm = forwardNet3(in_Optim,netxy,netz,mu);
    out_Optim_ML_norm = forwardNet3_xyzEnhance(in_Optim,netxy,netz,mu);
end
out_Optim_ML = out_Optim_ML_norm.*cat(3,sigx,sigy,sigz)+cat(3,mux,muy,muz);

% optim_FE
if ~isempty(csv_Optim)
    out_Optim_FE = readFEA_convBC_frame(csv_Optim,5,3); % frame=5, num_per_voxe=3.
else
    out_Optim_FE = [];
end

% whether use oriBC
if strcmp(BCType,'oriBC')
    outIntuiOri = oriBC(outIntuiOri);
    out_Optim_ML = oriBC(out_Optim_ML);
    if ~isempty(csv_Optim)
        out_Optim_FE = oriBC(out_Optim_FE);
    end
end

% ====== Shape comparison plots ======
% only type 3 stored here
colors = tab10;
surf(outIntuiOri(:,:,1),outIntuiOri(:,:,2),outIntuiOri(:,:,3), ...
    'FaceColor',0.99*[1,1,1],'FaceAlpha',0.8,'EdgeColor','none'); hold on;
error1 = 100*errorCalculate(outIntuiOri,out_Optim_ML,'relative1');
surf(out_Optim_ML(:,:,1),out_Optim_ML(:,:,2),out_Optim_ML(:,:,3),error1,...
    'FaceColor','none','FaceAlpha',0.8, 'EdgeColor',colors(1,:),'LineWidth',1.5); hold on;
% cb=colorbar; cb.Label.String='Error(%)';
clim([0, max(error1(:))])
if ~isempty(csv_Optim)
    surf(out_Optim_FE(:,:,1),out_Optim_FE(:,:,2),out_Optim_FE(:,:,3), ...
        'FaceColor','none','FaceAlpha',0.8,'EdgeColor',colors(2,:),'LineWidth',1.5); hold on;
end
axis equal; axis off;
light; camlight; material([0.4 0.5 0.7]); lighting gouraud;
colormap(viridis)
set(gca,'FontName','Arial','FontSize',18);

% Error plots
errorMap(outIntuiOri,out_Optim_ML,'relative1'); axis off;
if ~isempty(csv_Optim)
    errorMap(outIntuiOri,out_Optim_FE,'relative1'); axis off;
end

end

%%
function xyzPred = forwardNet3(ga_x,netD2Sxy,netD2Sz,mu)
ga_ins = reshape(ga_x,[15,15,2]);
ga_ins = (ga_ins-mu)/0.5; % normalization for ins (mu is not exactly 0.5)
xyPred = predict(netD2Sxy,ga_ins,'MiniBatchSize',1);
zPred = predict(netD2Sz,ga_ins,'MiniBatchSize',1);
xyzPred = cat(3,xyPred,zPred);
end

function dlin_round = round_of_dlin(dlin)
dlin_opti = tanh(dlin);
dlin_round = dlin_opti;
dlin_round(dlin_round>0) = 1; dlin_round(dlin_round<0) = -1; 
dlin_round = extractdata(dlin_round);
if isgpuarray(dlin_round)
    dlin_round = gather(dlin_round);
end
end

function xyzPred = forwardNet3_xyzEnhance(ga_x,netD2Sxy,netD2Sz,mu)
ga_ins = reshape(ga_x,[15,15,2]);
ga_ins = (ga_ins-mu)/0.5; % normalization for ins (mu is not exactly 0.5)
xyPred = predict(netD2Sxy,ga_ins,'MiniBatchSize',1);
zPred = predict(netD2Sz,ga_ins,'MiniBatchSize',1);

ga_ins_flipZ = flip(ga_ins,3);
ga_ins_temp2 = permute(ga_ins, [2,1,3,4]);
ga_ins_temp3 = flip(ga_ins_temp2,3);

xyPred_flipZ = predict(netD2Sxy,ga_ins_flipZ);
xyPred_temp2 = flip( permute(predict(netD2Sxy,ga_ins_temp2),[2,1,3,4]), 3);
xyPred_temp3 = flip( permute(predict(netD2Sxy,ga_ins_temp3),[2,1,3,4]), 3);
xyPred_mean = (xyPred + xyPred_flipZ + xyPred_temp2 + xyPred_temp3)/4;

zPred_flipZ = -predict(netD2Sz,ga_ins_flipZ);
zPred_temp2 = permute( predict(netD2Sz,ga_ins_temp2),[2,1,3,4] );
zPred_temp3 = permute( -predict(netD2Sz,ga_ins_temp3),[2,1,3,4] );
zPred_mean = (zPred + zPred_flipZ + zPred_temp2 + zPred_temp3)/4;

xyzPred = cat(3,xyPred_mean,zPred_mean);
end

function shape_oriBC = oriBC(shape)
[tempx,tempy,tempz,~] = data_convertBC1_back(shape(:,:,1),shape(:,:,2),shape(:,:,3),1);
shape_oriBC = cat(3,tempx,tempy,tempz);
end

% Compare shapes by R2 (R_squared) values
function R2 = getR2(xtrue,xpred)
% R2=getR2(xtrue,xpred)
temp1z = xtrue(:);
temp2z = xpred(:);
SSres = sum((temp2z-temp1z).^2);
SStot = sum((temp2z-mean(temp2z)).^2);
R2 = 1-SSres/SStot;
% R2 = 1-SSres/SStot*(length(temp1z)-1)/length(temp1z); % adjusted
end



