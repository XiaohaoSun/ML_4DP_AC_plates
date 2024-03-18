function xyzPred = forwardNet3(ga_x,netD2Sxy,netD2Sz,mu)
ga_ins = reshape(ga_x,[15,15,2]);
ga_ins = (ga_ins-mu)/0.5; % normalization for ins (mu is not exactly 0.5)
xyPred = predict(netD2Sxy,ga_ins,'MiniBatchSize',1);
zPred = predict(netD2Sz,ga_ins,'MiniBatchSize',1);
xyzPred = cat(3,xyPred,zPred);
end