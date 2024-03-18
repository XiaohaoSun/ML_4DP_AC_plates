function visualizeDesign(ins)


if length(size(ins)) == 3
    ins = [ins(:,:,1); ins(:,:,2)];
elseif size(ins,1) == 1
    ins = reshape(ins,15,15,2);
    ins = [ins(:,:,1); ins(:,:,2)];
end


%{
reverse in y-axis (1st dimension), 
so as to be compatible with surface plot and error map.
%}

% % 3D, to modify later
% voxImg = cat(3, ins(15:-1:1,:),ins(end:-1:16,:));
% 
% pprop.edgecolor             = 'w';
% pprop.backfacelighting      = 'unlit'; % 'reverselit' | 'unlit' | 'lit'.
% pprop.DiffuseStrength       = 0.9;
% pprop.ambientstrength       = 0.5;
% pprop.SpecularStrength      = 1;
% pprop.SpecularExponent      = 1;
% pprop.FaceLighting          = 'gouraud';
% 
% VOXview(voxImg, ones(size(voxImg))*0.9, 'colormap',parula, 'patch_props',pprop);
% % lpz = [-0, -1, 0.5]*size(voxImg, 1);
% % light('Position', lpz,'Style','local');
% 
% setwhite(gca)
% axis tight
% axis off


% === 2D ===
% colors = [0.07, 0.62, 1.00; 1.00, 0.41, 0.16; 0.39, 0.83, 0.07; 0.72, 0.27, 1.00]; % 蓝橙绿紫 
colors = [0.07, 0.62, 1.00; 1.00, 0.41, 0.16];

figure();
% top
subplot(2,1,1); imshow(ins(end:-1:16,:)); colormap(colors);
% bottom
subplot(2,1,2); imshow(ins(15:-1:1,:)); colormap(colors);

end