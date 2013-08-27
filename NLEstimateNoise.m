function [var] = NLEstimateNoise( image, winrad, nlrad )

if( ischar(image) )
    I = double(imread(image))./255;
else
    I = image;
end

[rows, cols, chans] = size(I);
imlen = rows*cols;
Igray = mean(I,3);
Igray = reshape( Igray, imlen, 1 );

% Assume 0 < var <= 1.
% So,  -Inf < varExp <= 0.
lastVarExp = 0;
varExp = 0;
A = NLAdjacency( image, winrad, nlrad, 1 );
for i=1:5
    var = 2^(varExp);
    
    % Adjust the matrix to reflect different denoising level.
    B = A.^(1/var);
    % Create Laplacian.
    L = speye(size(B)) + B;
    L = spdiags( sum(L,2), 0, imlen, imlen) \ L;
    % Noise matrix.
    N = speye(size(B)) - L;
    
    noise = reshape( full(N*Igray), rows, cols );
    % Trim to the actual denoised portion.
    noise = noise( (1+winrad):(rows-winrad), (1+winrad):(cols-winrad) );
    noise_spec = dct2(noise);
    % Sparsity is the wrong measure. Need to find
    % flatness of the spectrum.
    sparsity = sum(sum(abs(noise_spec))) / (size(noise,1)*size(noise,2));
    
    fprintf(1, 'VarExp: %d, Sparsity: %d\n', varExp, sparsity);
    varExp = varExp - 1;
end