function [ K_eqn, k1_eqn, k2_eqn ] = curvature_Analysis(Cx, Cy, Cz, inc, name, crop)

syms u v

[Nu, Nv] = size( Cx );
Nu = Nu - 1;
Nv = Nv - 1;

U = u.^(Nu:-1:0);

% V = v.^(Nv:-1:0);
VT = v.^( (Nv:-1:0)' );

Bu = BezierBasisMat( Nu );
Bv = BezierBasisMat( Nv );

Px = U*Bu*Cx*Bv'*VT;
Py = U*Bu*Cy*Bv'*VT;
Pz = U*Bu*Cz*Bv'*VT;

% Px = U*Bu*Cx*Bv'*V';
% Py = U*Bu*Cy*Bv'*V';
% Pz = U*Bu*Cz*Bv'*V';

sig = [Px; Py; Pz];

sig_u = diff( sig, u );
sig_v = diff( sig, v );

sig_uT = reshape( sig_u, 1, 3 );
sig_vT = reshape( sig_v, 1, 3 );

n = cross( sig_u, sig_v );
n = n ./ norm( n );

sig_uu = diff( sig_u, u );
sig_uv = diff( sig_u, v );
sig_vv = diff( sig_v, v );

sig_uuT = reshape( sig_uu, 1, 3 );
sig_uvT = reshape( sig_uv, 1, 3 );
sig_vvT = reshape( sig_vv, 1, 3 );

E = sig_uT * sig_u;
F = sig_uT * sig_v;
G = sig_vT * sig_v;

L = sig_uuT * n;
M = sig_uvT * n;
N = sig_vvT * n;

a = E * G - F^2;
b = -(L*G - 2*M*F + N*E);
c = (L * N - M^2);

K = c / a;

k1 = (-b - sqrt(b^2 - 4*a*c) ) / (2*a);
k2 = (-b + sqrt(b^2 - 4*a*c) ) / (2*a);

k1_eqn = matlabFunction( k1 );
k2_eqn = matlabFunction( k2 );

K_eqn = matlabFunction( K );

% dK_du = diff( K, u );
% dK_dv = diff( K, v );
% dK_dudv = diff( dK_du, v );

% dK_du_eqn = matlabFunction( dK_du );
% dK_dv_eqn = matlabFunction( dK_dv );
% dK_dudv_eqn = matlabFunction( dK_dudv );

if ~isempty(name)
    U_i = linspace(0, 1, inc)' * ones(1, inc);
    U_i = U_i( (1+crop):end-crop, (1+crop):end-crop );
    
    V_i = U_i';
    
    K_i = K_eqn(U_i, V_i);
    K_i( log10( abs( K_i ) ) < -8 ) = 0;
    close all
    figure()
    hold on
    s1 = surf( K_i );
    s1.EdgeColor = 'none';
    view( [-1, -1, 1] );
    title( strcat( name, ": Gaussian Curvature, K = \kappa_{1}\kappa_{2}" ) );
    hold off
    
    k1_i = k1_eqn(U_i, V_i);
    k1_i( log10( abs( k1_i ) ) < -8 ) = 0;
    figure()
    hold on
    s2 = surf( k1_i );
    s2.EdgeColor = 'none';
    view( [-1, -1, 1] );
    title( strcat( name, ": First Principal Curvature, \kappa_{1}") );
    hold off
    
    k2_i = k2_eqn(U_i, V_i);
    k2_i( log10( abs( k2_i ) ) < -8 ) = 0;
    figure()
    hold on
    s3 = surf( k2_i );
    s3.EdgeColor = 'none';
    view( [-1, -1, 1] );
    title( strcat( name, ": Second Principal Curvature, \kappa_{2}") );
    hold off
    
    figure()
    hold on
    s3 = surf( (k1_i + k2_i) / 2 );
    s3.EdgeColor = 'none';
    view( [-1, -1, 1] );
    title( strcat( name, ": Mean Curvature, (\kappa_{1} + \kappa_{2})/2") );
    hold off
    
    % dK_du_i = dK_du_eqn(U_i, V_i);
    % dK_du_i( log10( abs( dK_du_i ) ) < -8 ) = 0;
    % figure(4)
    % hold on
    % s4 = surf( dK_du_i );
    % s4.EdgeColor = 'none';
    % view( [-1, -1, 1] );
    % title( 'dK/du' );
    % hold off
    %
    % dK_dv_i = dK_dv_eqn(U_i, V_i);
    % dK_dv_i( log10( abs( dK_dv_i ) ) < -8 ) = 0;
    % figure(5)
    % hold on
    % s5 = surf( dK_dv_i );
    % s5.EdgeColor = 'none';
    % view( [-1, -1, 1] );
    % title( 'dK/dv' );
    % hold off
    %
    % dK_dudv_i = dK_dudv_eqn(U_i, V_i);
    % dK_dudv_i( log10( abs( dK_dudv_i ) ) < -8 ) = 0;
    % figure(6)
    % hold on
    % s6 = surf( dK_dudv_i );
    % s6.EdgeColor = 'none';
    % view( [-1, -1, 1] );
    % title( 'dK^{2}/dudv' );
    % hold off
end

end

function B = BezierBasisMat( n )

bi = @(i) factorial(n) / ( factorial(i) .* factorial(n - i) );

B = zeros( n+1 );

for i = 0:n
    b = bi(i) * Pascal_Level( n - i ); 
    ind = mod(n - i + 1, 2) + 1; 
    b(ind:2:end) = -1*b(ind:2:end); 
    B(i+1, 1:(n - i + 1) ) = b; 
end

end

function [pN] = Pascal_Level(N)

pN = 1;
for i = 1:N
    
    pN = [pN, 0] + [0, pN];

end

end