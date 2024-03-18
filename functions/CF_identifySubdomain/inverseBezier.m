function [Cx, Cy, Cz] = inverseBezier(U, V, Px, Py, Pz)

%{

U = [ ...; u_i^n, u_i^(n-1),..., u_i^2, u_i, 1;... ]
V = [ ...; v_j^m, v_j^(m-1),..., v_j^2, v_j, 1;... ]

length(U) == length(V), but m == n or m != n

%}

[~, Nu] = size(U);
Bu = BezierBasisMat( Nu - 1 );

[~, Nv] = size(V);
Bv = BezierBasisMat( Nv - 1 );

U_inv = pinv( U );
V_inv = pinv( V ) ;

Cx = ( Bu \ U_inv ) * Px * ( Bv \ V_inv )';
Cy = ( Bu \ U_inv ) * Py * ( Bv \ V_inv )';
Cz = ( Bu \ U_inv ) * Pz * ( Bv \ V_inv )';

end