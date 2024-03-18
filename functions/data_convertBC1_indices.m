function [outx,outy,outz] = data_convertBC1_indices(outx,outy,outz,indices)
% [outx,outy,outz] = data_convertBC1n(outx,outy,outz,indices), 
% indices -> [iy,ix], use iy=ix=2 by default.
% In the output data, say outx(:,:,1,1) 
% a column represents a y-direction material strip in FE model,
% a row represents a x-direction material strip in FE model.

if nargin == 3
    iy = 2;
    ix = 2;
else
    if length(indices) == 1
        iy = indices;
        ix = indices;
    else
        iy = indices(1);
        ix = indices(2);
    end
end

% indices for RP, i.e., reference point (0,0)
[iy0,~] = find(outy==0,1);
[~,ix0] = find(outx==0,1);

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
gamma = -atan2d(outy(iy,ix,:,:),outx(iy,ix,:,:))
temp1 = cosd(gamma).*outx - sind(gamma).*outy;
temp2 = sind(gamma).*outx + cosd(gamma).*outy;
outx = temp1;
outy = temp2;

% 2
gamma = atan2d(outz(iy,ix,:,:),outx(iy,ix,:,:))
temp1 = cosd(gamma).*outx + sind(gamma).*outz;
temp2 = -sind(gamma).*outx + cosd(gamma).*outz;
outx = temp1;
outz = temp2;

% 3
% now AC is x-axis, project BD on y-z plane to find the angle, then x-axis rotate
gamma = -atan2d(outz(iy,ix0,:,:)-outz(iy0,ix,:,:),outy(iy,ix0,:,:)-outy(iy0,ix,:,:))
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


end