function [outx,outy,outz,inBinary] = data_swapXY(outx,outy,outz,inBinary)

% In the output data, say outx(:,:,1,1) 
% a column represents a y-direction material strip in FE model,
% a row represents a x-direction material strip in FE model.

% 1
inBinary = permute(inBinary,[2,1,3,4]);

temp1 = permute(outx,[2,1,3,4]);
temp2 = permute(outy,[2,1,3,4]);
outx = temp2;
outy = temp1;
outz = permute(outz,[2,1,3,4]);

end
