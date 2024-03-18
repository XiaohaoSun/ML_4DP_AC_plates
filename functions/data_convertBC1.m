function [outx,outy,outz,inBinary] = data_convertBC1(outx,outy,outz,inBinary)

% In the output data, say outx(:,:,1,1) 
% a column represents a y-direction material strip in FE model,
% a row represents a x-direction material strip in FE model.

% if length(size(inBinary))==4 && length(size(outx))==4 % seems not necessary

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
gamma = -atan2d(outy(2,2,:,:),outx(2,2,:,:));
temp1 = cosd(gamma).*outx - sind(gamma).*outy;
temp2 = sind(gamma).*outx + cosd(gamma).*outy;
outx = temp1;
outy = temp2;

% 2
gamma = atan2d(outz(2,2,:,:),outx(2,2,:,:));
temp1 = cosd(gamma).*outx + sind(gamma).*outz;
temp2 = -sind(gamma).*outx + cosd(gamma).*outz;
outx = temp1;
outz = temp2;

% 3
% now AC is x-axis, project BD on y-z plane to find the angle, then x-axis rotate
gamma = -atan2d(outz(2,1,:,:)-outz(1,2,:,:),outy(2,1,:,:)-outy(1,2,:,:));
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

% % below is discarded
% % 1, z-axis rotate, AC to (1,1,0)-z plane
% % 2, (-1,1,0)-axis rotate, AC to (1,1,0)-axis
% % 3, (1,1,0)-axis rotate, B,D to same z level (angle not found yet, so not finished)
% gamma = 45-atan2d(outy(2,2,:,:),outx(2,2,:,:));
% temp1 = cosd(gamma).*outx - sind(gamma).*outy;
% temp2 = sind(gamma).*outx + cosd(gamma).*outy;
% outx = temp1;
% outy = temp2;
% 
% gama = atan2(outz(2,2,:,:),sqrt(2)*outx(2,2,:,:));
% rotm = axang2rotm([ones(length(gama(:)),1)*[-1,1,0],gama(:)]);
% for idx = 1:size(rotm,3)
%     temp1 = rotm(1,1,idx)*outx(:,:,:,idx) + rotm(1,2,idx)*outy(:,:,:,idx)...
%         + rotm(1,3,idx)*outz(:,:,:,idx);
%     temp2 = rotm(2,1,idx)*outx(:,:,:,idx) + rotm(2,2,idx)*outy(:,:,:,idx)...
%         + rotm(2,3,idx)*outz(:,:,:,idx);
%     temp3 = rotm(3,1,idx)*outx(:,:,:,idx) + rotm(3,2,idx)*outy(:,:,:,idx)...
%         + rotm(3,3,idx)*outz(:,:,:,idx);
%     outx(:,:,:,idx) = temp1;
%     outy(:,:,:,idx) = temp2;
%     outz(:,:,:,idx) = temp3;
% end


end