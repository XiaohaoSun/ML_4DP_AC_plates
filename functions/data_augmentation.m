function [outx,outy,outz,inBinary] = data_augmentation(outx,outy,outz,inBinary)

% In the output data, say outx(:,:,1,1) 
% a column represents a y-direction material strip in FE model,
% a row represents a x-direction material strip in FE model.

[temp1,temp2,temp3,temp4] = data_flipX(outx,outy,outz,inBinary);
outx = cat(4,outx,temp1);
outy = cat(4,outy,temp2);
outz = cat(4,outz,temp3);
inBinary = cat(4,inBinary,temp4);

[temp1,temp2,temp3,temp4] = data_flipY(outx,outy,outz,inBinary);
outx = cat(4,outx,temp1);
outy = cat(4,outy,temp2);
outz = cat(4,outz,temp3);
inBinary = cat(4,inBinary,temp4);

[temp1,temp2,temp3,temp4] = data_swapXY(outx,outy,outz,inBinary);
outx = cat(4,outx,temp1);
outy = cat(4,outy,temp2);
outz = cat(4,outz,temp3);
inBinary = cat(4,inBinary,temp4);

[temp1,temp2,temp3,temp4] = data_flipZ(outx,outy,outz,inBinary);
outx = cat(4,outx,temp1);
outy = cat(4,outy,temp2);
outz = cat(4,outz,temp3);
inBinary = cat(4,inBinary,temp4);

end
