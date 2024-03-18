function [outx,outy,outz,inBinary] = data_convertBCori(outx,outy,outz,inBinary)
% Copy of data_convertBC1_back.

% In the output data, say outx(:,:,1,1) 
% a column represents a y-direction material strip in FE model,
% a row represents a x-direction material strip in FE model.

% if length(size(inBinary))==4 && length(size(outx))==4

% Sequentially using x,y or z rotation enables vectorized operation,
% which is faster than using rotm by axang2rotm which uses loops.
% 1, z-axis rotate, AC to xz-plane
% 2, y-axis rotate, AC to x-axis
% 3, x-axis rotate, B,D to same z level
% 4, z-axis rotate, AC to (1,1,0)-axis

%******************
% *****B-----C*****
% *****|     |*****
% *****|     |*****
% *****A-----D*****
% *****************

% 1
gamma = -atan2d(outy(end,end,:,:),outx(end,end,:,:));
temp1 = cosd(gamma).*outx - sind(gamma).*outy;
temp2 = sind(gamma).*outx + cosd(gamma).*outy;
outx = temp1;
outy = temp2;

% 2
gamma = atan2d(outz(end,end,:,:),outx(end,end,:,:));
temp1 = cosd(gamma).*outx + sind(gamma).*outz;
temp2 = -sind(gamma).*outx + cosd(gamma).*outz;
outx = temp1;
outz = temp2;

% 3
% now AC is x-axis, project BD on y-z plane to find the angle, then x-axis rotate
gamma = -atan2d(outz(end,1,:,:)-outz(1,end,:,:),outy(end,1,:,:)-outy(1,end,:,:));
temp1 = cosd(gamma).*outy - sind(gamma).*outz;
temp2 = sind(gamma).*outy + cosd(gamma).*outz;
outy = temp1;
outz = temp2;

% 4
gamma = 45;
temp1 = cosd(gamma).*outx - sind(gamma).*outy;
temp2 = sind(gamma).*outx + cosd(gamma).*outy;
outx = temp1;
outy = temp2;
    
% end
    
end