function [Idenoised, noise] = NLDenoise( image, winrad, nlrad, simvar )

if( ischar(image) )
    I = double(imread(image))./255;
else
    I = image;
end

[rows, cols, chans] = size(I);
imlen = rows*cols;
Igray = mean(I,3);

L = NLWeights( NLAdjacency(image, winrad, nlrad, simvar) );

%Igray = reshape(Igray,imlen,1);
%inds = reshape(1:imlen,rows,cols);
I = reshape(I((1+winrad):(rows-winrad),(1+winrad):(cols-winrad),:), (rows-2*winrad)*(cols-2*winrad), chans);

Idenoised = reshape( full(L*I), (rows-2*winrad), (cols-2*winrad), chans );
noise = reshape(I,(rows-2*winrad),(cols-2*winrad), chans) - Idenoised;