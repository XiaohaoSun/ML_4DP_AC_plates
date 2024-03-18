function [n, m, Cx, Cy, Cz, err] = Best_Bezier(Px, Py, Pz, U_ind, V_ind, doPlot)

[Nu, Nv] = size(Px);
if isempty(U_ind) || isempty(V_ind) 
    U_ind = linspace(0, 1, Nu)';
    V_ind = linspace(0, 1, Nv)';
end

interp_x = zeros( Nu - 1, Nv - 1 );
interp_y = zeros( Nu - 1, Nv - 1 );
interp_z = zeros( Nu - 1, Nv - 1 );
for i = 1:Nu-1
    for j = 1:Nv-1
        interp_x(i, j) = sum( sum( Px( i:(i+1), j:(j+1) ) ) ) / 4;
        interp_y(i, j) = sum( sum( Py( i:(i+1), j:(j+1) ) ) ) / 4;
        interp_z(i, j) = sum( sum( Pz( i:(i+1), j:(j+1) ) ) ) / 4;
    end
end
U_interp = ( U_ind(1:end-1) + U_ind(2:end) ) / 2;
V_interp = ( V_ind(1:end-1) + V_ind(2:end) ) / 2;

u_range = [2, min( floor(Nu/3), 10) ];
v_range = [2, min( floor(Nv/3), 10) ];

n = 0;
m = 0;
Cx = 0;
Cy = 0;
Cz = 0;
err = 1e6;

for i = u_range(1):u_range(2)
    for j = v_range(1):v_range(2)
        p_n = @(u) u.^(i:-1:0);
        p_m = @(v) v.^(j:-1:0);
        U = p_n( U_ind );
        V = p_m( V_ind );
        
        [Cx_r, Cy_r, Cz_r] = inverseBezier(U, V, Px, Py, Pz);
        
        Bu = BezierBasisMat( i );
        Bv = BezierBasisMat( j );
        
        U_int = p_n( U_interp );
        V_int = p_m( V_interp );
        
        Px_int = U_int*Bu*Cx_r*Bv'*V_int';
        Py_int = U_int*Bu*Cy_r*Bv'*V_int';
        Pz_int = U_int*Bu*Cz_r*Bv'*V_int';

        E_int_x = abs( Px_int - interp_x );
        E_int_y = abs( Py_int - interp_y );
        E_int_z = abs( Pz_int - interp_z );

        E = sum(E_int_x(:)) + sum(E_int_y(:)) + sum(E_int_z(:));
%         Px_r = U*Bu*Cx_r*Bv'*V';
%         Py_r = U*Bu*Cy_r*Bv'*V';
%         Pz_r = U*Bu*Cz_r*Bv'*V';
%         E_x = abs( Px - Px_r );
%         E_y = abs( Py - Py_r );
%         E_z = abs( Pz - Pz_r );  
%         E = sum(E_x(:)) + sum(E_y(:)) + sum(E_z(:));
        if E < err
            err = E;
            n = i;
            m = j;
            Cx = Cx_r;
            Cy = Cy_r;
            Cz = Cz_r;
        end
    end
end

if Cx ~= 0 & doPlot ~= 0
    
    close all
    
    figure(1)
    hold on
    s1 = surf( Px, Py, Pz );
    s1.EdgeColor = 'none';
    view( [1, 1, 1] );
    axis equal
    hold off

    U = U_ind.^(n:-1:0);
    V = V_ind.^(m:-1:0);
    
    Bu = BezierBasisMat( n );
    Bv = BezierBasisMat( m );
    Px_b = U*Bu*Cx*Bv'*V';
    Py_b = U*Bu*Cy*Bv'*V';
    Pz_b = U*Bu*Cz*Bv'*V';
    
    figure(2)
    hold on
    s2 = surf( Px_b, Py_b, Pz_b );
    s2.EdgeColor = 'none';
    view( [1, 1, 1] );
    axis equal
    hold off
    
    figure(3)
    hold on
    s3 = surf( abs(Px_b - Px) + abs(Py_b - Py) + abs(Pz_b - Pz) );
    s3.EdgeColor = 'none';
    view( [1, 1, 1] );
    hold off
    
    figure(4)
    hold on
    s4 = surf( (abs(Px_b - Px) + abs(Py_b - Py) + abs(Pz_b - Pz))./...
                norm( Px + Py + Pz ) );
    s4.EdgeColor = 'none';
    view( [1, 1, 1] );
    hold off
    
%     figure(4)
%     hold on
%     s4 = surf( U_i * ones(1, Nv), ones(Nu, 1) * V_i', abs(Pz_b - Pz) );
%     s4.EdgeColor = 'none';
%     axis equal
%     view( [1, 1, 1] );
%     hold off
    
end


end




