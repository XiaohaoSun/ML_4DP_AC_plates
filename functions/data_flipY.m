function [outx,outy,outz,inBinary] = data_flipY(outx,outy,outz,inBinary)

% In the output data, say outx(:,:,1,1) 
% a column represents a y-direction material strip in FE model,
% a row represents a x-direction material strip in FE model.


% 1, data order rearrangement
% 2, shape transformation

% 1
inBinary = inBinary(end:-1:1,:,:,:);
outx = outx(end:-1:1,:,:,:);
outy = outy(end:-1:1,:,:,:);
outz = outz(end:-1:1,:,:,:);

% 2
outx = outx - outx(1,1,:,:);
outy = outy - outy(1,1,:,:);
outz = outz - outz(1,1,:,:);
outy = -outy;

gamma = 45-atan2d(outy(end,end,:,:),outx(end,end,:,:));
temp1 = cosd(gamma).*outx - sind(gamma).*outy;
temp2 = sind(gamma).*outx + cosd(gamma).*outy;
outx = temp1;
outy = temp2;

end
