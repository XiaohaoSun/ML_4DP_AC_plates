function outXYZ = readFEA_convBC_frame(csvFile,frame,num_per_voxel)

dataFE = readmatrix(csvFile);
dataFE = sortrows(dataFE,[3 2]);

c_x_2D = reshape(dataFE(:,2)+dataFE(:,4+(frame-1)*3+1),[46,46])';
c_y_2D = reshape(dataFE(:,3)+dataFE(:,4+(frame-1)*3+2),[46,46])';
c_z_2D = reshape(dataFE(:,4+(frame-1)*3+3),[46,46])';

inBinary = 1; 
[outx,outy,outz,inBinary] = data_convertBC1(...
    c_x_2D,c_y_2D,c_z_2D,inBinary);

outXYZ = cat(3,outx,outy,outz);

if nargin==3
    outXYZ = outXYZ((1):num_per_voxel:end,...
        (1):num_per_voxel:end,:); % choose corner-point
end

end