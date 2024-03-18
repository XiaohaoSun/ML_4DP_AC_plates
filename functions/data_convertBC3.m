function [outx,outy,outz,inBinary] = data_convertBC3(outx,outy,outz,inBinary)

% In the output data, say outx(:,:,1,1) 
% a column represents a y-direction material strip in FE model,
% a row represents a x-direction material strip in FE model.

% Sequentially using x,y or z rotation enables vectorized operation,
% which is faster than using rotm by axang2rotm which uses loops.
% 1, z-axis rotate, AB to xz-plane
% 2, y-axis rotate, AB to x-axis
% 3, x-axis rotate, C to xy-plane(i.e., z=0, same level with B)

%******************
% *****D-----C*****
% *****|     |*****
% *****|     |*****
% *****A-----B*****
% *****************

% 1
gamma = -atan2d(outy(1,end,:,:),outx(1,end,:,:));
temp1 = cosd(gamma).*outx - sind(gamma).*outy;
temp2 = sind(gamma).*outx + cosd(gamma).*outy;
outx = temp1;
outy = temp2;

% 2
gamma = atan2d(outz(1,end,:,:),outx(1,end,:,:));
temp1 = cosd(gamma).*outx + sind(gamma).*outz;
temp2 = -sind(gamma).*outx + cosd(gamma).*outz;
outx = temp1;
outz = temp2;

% 3
% now AB is x-axis, project BC on yz-plane to find the angle, then x-axis rotate
gamma = -atan2d(outz(end,end,:,:),outy(end,end,:,:));
temp1 = cosd(gamma).*outy - sind(gamma).*outz;
temp2 = sind(gamma).*outy + cosd(gamma).*outz;
outy = temp1;
outz = temp2;

    
end