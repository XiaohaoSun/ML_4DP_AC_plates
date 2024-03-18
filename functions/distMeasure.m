function [distTable, distML] = distMeasure(target, optim, refL)
% target can be triangulation or 15x15x3 array, optim is the achieved shape.

if nargin == 2
    refL = 40; % default value is 40
end

if ~strcmp(class(target),'triangulation') % not triangulation
    TOx = double(target(:,:,1));
    TOy = double(target(:,:,2));
    TOz = double(target(:,:,3));
    % TOT = delaunay(TOx,TOy);
    TOT = getConnectivity(TOx,TOy,TOz);
    target = triangulation(TOT,TOx(:),TOy(:),TOz(:));
end

QP1 = double(reshape(optim,[],3));
distML = point2trimesh('Faces',target.ConnectivityList, ...
    'Vertices',target.Points, 'QueryPoints',QP1, 'Algorithm','vectorized');

% dist map
figure();
distML_perc = abs(distML)/refL*100;
s = pcolor(reshape(distML_perc,16,16));
cb=colorbar; 
cb.Label.String='Relative distance(%)';
s.FaceColor = 'interp';
s.EdgeColor = 'w';
title(['Max ',num2str(max(distML_perc),3),'%, Avg ',num2str(mean(distML_perc),3),'%']);
set(gca,'FontName','Arial','FontSize',15); axis off;
colormap(viridis);

% dist table  
Objects = ["Target vs. Optim"];
distAvg(1) = mean(abs(distML));
distMax(1) = max(abs(distML));
distRMS(1) = rms(distML);
distTable = table(Objects,distAvg,distMax,distRMS);

% Final reshape
distML = reshape(distML,16,16);

end