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