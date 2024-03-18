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