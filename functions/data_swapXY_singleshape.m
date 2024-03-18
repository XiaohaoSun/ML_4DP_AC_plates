function shapeOut = data_swapXY_singleshape(shapeIn)
% chirality will be changed.

% In the output data, say outx(:,:,1,1) 
% a column represents a y-direction material strip in FE model,
% a row represents a x-direction material strip in FE model.

outx = shapeIn(:,:,1);
outy = shapeIn(:,:,2);
outz = shapeIn(:,:,3);


temp1 = permute(outx,[2,1]);
temp2 = permute(outy,[2,1]);
outx = temp2;
outy = temp1;
outz = permute(outz,[2,1]);

shapeOut = cat(3,outx,outy,outz);

end
