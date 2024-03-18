function [mu,sig,out] = NormData(data)
% 4D data used for CNN
mu = mean(data(:));
sig =  std(data(:));
out = (data - mu)/sig;
end