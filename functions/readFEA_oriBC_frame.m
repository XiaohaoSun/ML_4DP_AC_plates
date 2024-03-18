function outXYZ = readFEA_oriBC_frame(csvFile,frame,num_per_voxel)
% csvFile is file name.
% frame can be simply set to 5 in current case. (frame 5 corresponds to the deformation at strain 0.05)
% num_per_voxel=1 gives 46x46x3 coordinates; =3 gives 16x16x3 coordinates.

dataFE = readmatrix(csvFile);
dataFE = sortrows(dataFE,[3 2]);

c_x_2D = reshape(dataFE(:,2)+dataFE(:,4+(frame-1)*3+1),[46,46])';
c_y_2D = reshape(dataFE(:,3)+dataFE(:,4+(frame-1)*3+2),[46,46])';
c_z_2D = reshape(dataFE(:,4+(frame-1)*3+3),[46,46])';

outXYZ = cat(3,c_x_2D,c_y_2D,c_z_2D);

if nargin==3
    outXYZ = outXYZ((1):num_per_voxel:end,...
        (1):num_per_voxel:end,:); % choose corner-point
end

end