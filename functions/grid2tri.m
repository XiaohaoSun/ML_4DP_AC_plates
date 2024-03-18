function target = grid2tri(target)

TOx = double(target(:,:,1));
TOy = double(target(:,:,2));
TOz = double(target(:,:,3));
% TOT = delaunay(TOx,TOy);
TOT = getConnectivity(TOx,TOy,TOz);
target = triangulation(TOT,TOx(:),TOy(:),TOz(:));

end