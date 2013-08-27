function [ constrained, vals ] = scribData( scribImage, winrad )
% Converts scribble image into usable constraints
%  constrained - constrained(i)==1 if pixel i is constrained. 0 otherwise.
%  vals        - vals(i)==0 if pixel i is unconstrained or constrained to 0. 1 otherwise.
%  scribimage  - Normalizes [0,1] image of scribbles. 0 if constrained to 0, 1 if constrained to 1.
%  winrad      - Ignore pixels closer than this to the edge

if( ~exist('winrad', 'var') || isempty(winrad) )
    winrad = 0;
end

if( ischar(scribImage) )
    scribs = double(imread(scribImage)) ./ 255;
else
    scribs = scribImage;
end

[rows,cols,chans] = size(scribs);
scribs = scribs((1+winrad):(rows-winrad),(1+winrad):(cols-winrad),1);

constrained = (scribs == 0);
constrained = constrained + (scribs == 1);

vals = constrained .* scribs;
