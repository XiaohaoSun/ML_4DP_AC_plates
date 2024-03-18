function [subspace] = find_deviation(ML, target, tol, flagCurvature)
% tol >= 1, subspace will have show the largest tol-number nodes
% 0 <= tol < 1, subspace will show all nodes with comparison > tol

if nargin==3
    flagCurvature = 0;
end

if flagCurvature == 1
    % 
    [r, c, ~] = size(ML);
    
    [~, ~, Cx_ML, Cy_ML, Cz_ML, ~] =...
          Best_Bezier(ML(:, :, 1), ML(:, :, 2), ML(:, :, 3), [], [], 0);
    
    [~, ~, Cx_T, Cy_T, Cz_T, ~] =...
          Best_Bezier(target(:, :, 1), target(:, :, 2), target(:, :, 3), [], [], 0);
      
    [ ~, k1_eqn_ML, k2_eqn_ML ] = curvature_Analysis(Cx_ML, Cy_ML, Cz_ML, 0, [], 0);
    [ ~, k1_eqn_T, k2_eqn_T ] = curvature_Analysis(Cx_T, Cy_T, Cz_T, 0, [], 0);
    
    inc = 16;
    U_i = linspace(0, 1, inc)' * ones(1, inc);
    V_i = U_i';
    
    k1_ML = k1_eqn_ML(U_i, V_i);
    k2_ML = k2_eqn_ML(U_i, V_i);
    
    k1_T = k1_eqn_T(U_i, V_i);
    k2_T = k2_eqn_T(U_i, V_i);
    
    % === Function Weights (Probably do not want to change theses
    % t = 'Dual-Sigmoid-esque';
    % A = 1; 
    % B1 = 8; B2 = B1;
    % C1 = 0.8; C2 = C1;
    % f = @(u, v) A ./ ( (1 + exp(B1.*(u - C1) ) ) .* (1 + exp(B2.*(v - C2) ) ) );
    
    % t = 'Parabolic-esque'
    % n = 0.75;
    % a = 1 / (0.5)^n / (1 - 0.5)^n;
    % f = @(u, v) a .* (u.^n .* (1 - u).^n) .* (v.^n .* (1 - v).^n);
    
    % t = 'Exclude Edges';
    % f = @(u, v) ceil( mod(u, 1) .* mod(v, 1) );
    
    t = 'No Exclusions';
    f = @(u, v) ones( size(u) );
    
    % === Comparison functions, either consider curvature or distance
    F = sqrt( (k1_ML - k1_T).^2 + (k2_ML - k2_T).^2 ) .* f(U_i, V_i);

elseif flagCurvature == 0
    % === Comparison functions, either consider curvature or distance
    F = sqrt( (ML(:,:,1)-target(:,:,1)).^2 +...
              (ML(:,:,2)-target(:,:,2)).^2 +...
              (ML(:,:,3)-target(:,:,3)).^2);
end

% F is either consider curvature or distance
compared = F / max( F(:) );

if tol < 1
    subspace = compared > tol;
else
    [~, ind] = sort( compared(:), 'descend' );
    subspace = zeros( size( compared ) );
    subspace( ind(1:tol) ) = 1;
end


end