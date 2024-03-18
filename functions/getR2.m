function R2 = getR2(trueValue, predValue)
% R2 = getR2( trueValue, predValue)

trueValue = trueValue(:);
predValue = predValue(:);

SSres = sum((predValue-trueValue).^2);
SStot = sum((trueValue-mean(trueValue)).^2);
R2 = 1-SSres/SStot;

end