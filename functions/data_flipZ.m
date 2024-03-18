function [outx,outy,outz,inBinary] = data_flipZ(outx,outy,outz,inBinary)

% In the output data, say outx(:,:,1,1) 
% a column represents a y-direction material strip in FE model,
% a row represents a x-direction material strip in FE model.

inBinary = inBinary(:,:,end:-1:1,:);
outz = -outz;

end
