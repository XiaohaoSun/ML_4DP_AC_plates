function visualize_TargetTri(targetIn, flag)
% flag=0, do nothting.
% flag=1, sample targetIn (46x46x3) to get a 16x16x3 shape.

if nargin == 1
    flag = 0;
end

if flag == 1
    targetIn = targetIn(1:3:end,1:3:end,:);
end

shape_Target = targetIn;

shape_Target_TR = grid2tri(shape_Target);

colors = tab20;

trisurf(shape_Target_TR,'FaceColor',colors(16,:),'FaceAlpha',0.5, ...
    'EdgeColor',colors(9,:),'LineWidth',1.5); 
% trisurf(shape_Target_TR,'FaceColor','none','EdgeColor',colors(9,:),'LineWidth',1.5); 

% Plot settings
axis equal; axis off;
% light('Position',[0,0,1],'Style','infinite');
camlight; 
% material([0.4 0.5 0.7]); 
% lighting gouraud;

end