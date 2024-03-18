function errorTable_return = errorTable(optimType,BCType, ...
    outIntui, in_Optim, csv_Optim, insIntui, ...
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
    out_Optim_ML_norm = forwardNet3(in_Optim,netxy,netz,mu);
end
out_Optim_ML = out_Optim_ML_norm.*cat(3,sigx,sigy,sigz)+cat(3,mux,muy,muz);

% optim_FE
if ~isempty(csv_Optim)
    out_Optim_FE = readFEA_convBC_frame(csv_Optim,5,3); % frame=5, num_per_voxe=3.
else
    out_Optim_FE = [];
end

% ML of insIntui
if ~isempty(insIntui)
    out_ML_insIntui_norm = forwardNet3(...
        cat(3,insIntui(1:15,:),insIntui(16:30,:)),netxy,netz,mu);
    out_ML_insIntui = out_ML_insIntui_norm.*cat(3,sigx,sigy,sigz)...
        +cat(3,mux,muy,muz);
end

% whether use oriBC
if strcmp(BCType,'oriBC')
    outIntuiOri = oriBC(outIntuiOri);
    out_Optim_ML = oriBC(out_Optim_ML);
    if ~isempty(csv_Optim)
        out_Optim_FE = oriBC(out_Optim_FE);
    end
    if ~isempty(insIntui)
        out_ML_insIntui = oriBC(out_ML_insIntui);
    end
end

% Error table for statistical quantification 
Objects = ["Target vs. ML-optim"; "Target vs. FE-optim"; "FE vs. ML, optim"; "FE vs. ML, intui"];
ErrAvg = NaN(4,1);
RMSE = NaN(4,1);
RMSEz = NaN(4,1);
R2 = NaN(4,1);
R2z = NaN(4,1);

% outIntuiOri, Target shape. Or FE of insIntui for intuitive design.
% out_Optim_ML, ML of optimal shape.
% out_Optim_FE, FE of optimal shape.
% out_ML_insIntui, ML of intuitive.

ErrAvg(1) = mean(errorCalculate(outIntuiOri, out_Optim_ML, 'absolute'),'all');
RMSE(1) = rms(outIntuiOri-out_Optim_ML,'all');
RMSEz(1) = rms(outIntuiOri(:,:,3)-out_Optim_ML(:,:,3),'all');
R2(1) = getR2(outIntuiOri,out_Optim_ML);
R2z(1) = getR2(outIntuiOri(:,:,3),out_Optim_ML(:,:,3));

if ~isempty(csv_Optim)
    ErrAvg(2) = mean(errorCalculate(outIntuiOri, out_Optim_FE, 'absolute'),'all');
    RMSE(2) = rms(outIntuiOri-out_Optim_FE,'all');
    RMSEz(2) = rms(outIntuiOri(:,:,3)-out_Optim_FE(:,:,3),'all');
    R2(2) = getR2(outIntuiOri,out_Optim_FE);
    R2z(2) = getR2(outIntuiOri(:,:,3),out_Optim_FE(:,:,3));

    ErrAvg(3) = mean(errorCalculate(out_Optim_FE, out_Optim_ML, 'absolute'),'all');
    RMSE(3) = rms(out_Optim_FE-out_Optim_ML,'all');
    RMSEz(3) = rms(out_Optim_FE(:,:,3)-out_Optim_ML(:,:,3),'all');
    R2(3) = getR2(out_Optim_FE,out_Optim_ML);
    R2z(3) = getR2(out_Optim_FE(:,:,3),out_Optim_ML(:,:,3));
end
if ~isempty(insIntui)
    % in this case, FE of insIntui is the target, i.e., outIntuiOri.
    ErrAvg(4) = mean(errorCalculate(outIntuiOri, out_ML_insIntui, 'absolute'),'all');
    RMSE(4) = rms(outIntuiOri-out_ML_insIntui,'all');
    RMSEz(4) = rms(outIntuiOri(:,:,3)-out_ML_insIntui(:,:,3),'all');
    R2(4) = getR2(outIntuiOri,out_ML_insIntui);
    R2z(4) = getR2(outIntuiOri(:,:,3),out_ML_insIntui(:,:,3));
end

errorTable_return = table(Objects,ErrAvg,RMSE,R2,RMSEz,R2z);

end

%%
function xyzPred = forwardNet3(ga_x,netD2Sxy,netD2Sz,mu)
ga_ins = reshape(ga_x,[15,15,2]);
ga_ins = (ga_ins-mu)/0.5; % normalization for ins (mu is not exactly 0.5)
xyPred = predict(netD2Sxy,ga_ins,'ExecutionEnvironment','gpu','MiniBatchSize',1);
zPred = predict(netD2Sz,ga_ins,'ExecutionEnvironment','gpu','MiniBatchSize',1);
xyzPred = cat(3,xyPred,zPred);
end

function dlin_round = round_of_dlin(dlin)
dlin_opti = tanh(dlin);
dlin_round = dlin_opti;
dlin_round(dlin_round>0) = 1; dlin_round(dlin_round<0) = -1; 
dlin_round = extractdata(dlin_round);
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
