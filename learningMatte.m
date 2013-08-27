function [ alpha ] = learningMatte( image, scribs )
%LEARNINGMATTE Makes the Laplacian given by Y. Zheng, C. Kambhamettu. "Learning Based Digital Matting". ICCV, 2009.
%   image  - Path to image
%   scribs - Path to scribbles

%A = NLAdjacency(image);

I = double(imread(image))./255;
scribs = double(imread(scribs)) ./255;
scribs = mean(scribs,3);

[rows, cols, chans] = size(I);
imsize = rows*cols;

[constr, cvals] = scribData( scribs );
constr = reshape( constr, imsize, 1 );
cvals = reshape( cvals, imsize, 1 );

inds = reshape([1:imsize], rows, cols);
I = reshape( I, imsize, chans );

% Some constants.
winRad = 2;
hoodsize = (2*winRad+1)^2;
lambda = 1e-2;
c = 1e9;

rowinds = zeros(rows*cols - (rows-2*winRad)*(cols-2*winRad),1);
hood = zeros( size(rowinds,1), hoodsize );
f = ones( size(rowinds,1), hoodsize);

k = 1;
for i=1+winRad:rows-winRad
    for j=1+winRad:cols-winRad
        hood( k, : ) = reshape( inds( i-winRad:i+winRad, j-winRad:j+winRad ), 1, hoodsize);
        rowinds( k, 1 ) = inds(i,j);
        
        X = I( hood(k,:), : );
        X(:, end+1) = ones(hoodsize,1);

        xprime = I( inds(i,j), : );
        xprime( :, end+1) = 1;
        f(k,:) = ((X*X' + lambda*eye(hoodsize) ) \ X) * xprime';
        
        k = k + 1;
    end
end

%F = sparse( repmat(1:hoodSize, imsize, 1), hood, f, imsize, imsize );
F = sparse( repmat( rowinds, 1, hoodsize), hood, f, imsize, imsize );
C = spdiags( c * constr, 0, imsize, imsize );
L2 = speye(imsize) - F;

alpha = ((L2'*L2 + C) \ C) * cvals;
alpha = reshape( alpha, rows, cols );
