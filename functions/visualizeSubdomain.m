function visualizeSubdomain(subspace)

%{
reverse in y-axis (1st dimension), 
so as to be compatible with surface plot and error map.
%}

% === 2D ===
% colors = [0.07, 0.62, 1.00; 1.00, 0.41, 0.16; 0.39, 0.83, 0.07; 0.72, 0.27, 1.00]; % 蓝橙绿紫 

colors = [0.44*ones(1,3); [255,222,23]/255];
figure(); imshow(sub2full(subspace(end:-1:1,:),[size(subspace)*20,1])); colormap(colors);

end