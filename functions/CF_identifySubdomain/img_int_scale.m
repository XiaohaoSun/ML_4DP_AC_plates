function img_scaled = img_int_scale(img, scale)

scale = round( scale );

[m, n, p] = size( img );
img_scaled = zeros( m*scale, n*scale, p );
pix = ones(scale, scale, p);
if isa( img(1,1,1), 'uint8' )
    pix = uint8( pix );
    img_scaled = uint8( img_scaled );
end
for i = 1:m
    for j = 1:n
        if sum( img(i, j, :) ) > 0
            I = ((i-1)*scale + 1):i*scale;
            J = ((j-1)*scale + 1):j*scale;
            img_scaled( I, J, : ) = img(i, j, :) .* pix;
        end
    end
end

end
