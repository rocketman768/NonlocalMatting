function [A, patches, patchNdx] = NLAdjacency( image, winrad, nlrad, simvar, maxNeighbors, doGray, doGrad, doPatches, distanceVar )
% Calculates the NL adjacency matrix.
% --Return values--
% A - The adjacency matrix. All entries in [0,1].
% patches - A (patchsize*channels)x(#patches) matrix containing
%           all the patches.
% patchNdx - 2x(#patches) matrix where each row corresponds to (x,y)
%            location of patch in image.
% --Arguments--
% image - either image name or image.
% winrad - The radius of the windows. Default=3.
% nlrad - The radius of the NL search window. Default=10.
% simvar - The variance of the similarity kernel. Default=2e-1.
% maxNeighbors - The maximum number of neighbors for each pixel. I.e., the
%                maximum # of nonzero entries in each row.
% doPatches - 1 if you want patches. 0 otherwise. Default=0.
% distanceVar - Variance of pixel distance term. Default=Inf.

if( ischar(image) )
    I = double(imread(image))./255;
else
    I = image;
end

if( ~exist('winrad') || isempty(winrad) )
    winrad = 3;
end
winsize = (2*winrad+1)^2;
if( ~exist('nlrad') || isempty(nlrad) )
    nlrad = 10;
end
nlsize = (2*nlrad+1)^2;
if( ~exist('simvar') || isempty(simvar) )
    simvar = 2e-1;
end
if( ~exist('doPatches') || isempty(doPatches) )
    doPatches = 0;
end
if( ~exist('maxNeighbors') || isempty(maxNeighbors) )
    maxNeighbors = nlsize;
end
if( ~exist('distanceVar') || isempty(distanceVar) )
    distanceVar = Inf;
end

if( ~exist('doGray') || isempty(doGray) )
    doGray = 0;
end
if( doGray == 1 )
    I = mean(I,3);
end

if( ~exist('doGrad') || isempty(doGrad) )
    doGrad = 0;
end
    xGrad = conv2(mean(I,3), [-.5,.5], 'same');
    gradVar = 1e0;

[rows, cols, chans] = size(I);
imlen = rows*cols;
centlen = (rows-2*winrad)*(cols-2*winrad); % Number of central pixels.

%blah = sqrt(2*winsize-1);
%sigma = sqrt(simvar);
h = fspecial('gaussian', [2*winrad+1,2*winrad+1], winrad);
h = reshape(h,winsize,1);

%inds = reshape([1:imlen], rows, cols);
inds = reshape( (1:centlen), rows-2*winrad, cols-2*winrad);
c = 1;
% s = -1*ones(1, ceil(nlsize/2)*imlen);
% rowinds = -1*ones(1, ceil(nlsize/2)*imlen);
% colinds = -1*ones(1, ceil(nlsize/2)*imlen);

%s = sparse(imlen,nlsize);
%rowinds = sparse(imlen,nlsize);
%colinds = sparse(imlen,nlsize);

s = -1*ones(1, maxNeighbors*imlen);
rowinds = -1*ones(1, maxNeighbors*centlen);
colinds = -1*ones(1, maxNeighbors*centlen);

if( doPatches )
    patches = zeros(winsize*chans, imlen);
    patchNdx = zeros(2,imlen);
end

for i=1+winrad:rows-winrad % Corresponds to original image pixel row.
    for j=1+winrad:cols-winrad % Corresponds to original image pixel column.
        % Find the adjacency weights to pixel (i,j)

        w = inds(i-winrad,j-winrad);
        %win_w_inds = reshape( inds( i-winrad:i+winrad, j-winrad:j+winrad ), winsize , 1);
        %win_w = Igray(win_w_inds);
        win_w = reshape( I( i-winrad:i+winrad, j-winrad:j+winrad, : ), winsize, chans);
        win_wGrad = reshape( xGrad( i-winrad:i+winrad, j-winrad:j+winrad, : ), winsize, 1 );

        if( doPatches )
            patches(:,c) = reshape(win_w, winsize*chans, 1);
            patchNdx(:,c) = [i;j];
        end

        %wins = zeros(winsize,1);
        %d = c;
        similarities = zeros(1,nlsize);
        simNdx = zeros(1,nlsize);
        pp = 1;
        for a= max(i-nlrad, 1+winrad):min(i+nlrad, rows-winrad)
            for b= max(j-nlrad, 1+winrad):min(j+nlrad, cols-winrad)

                %if( a == i && b == j )
                %    continue;
                %end

                k = inds(a-winrad,b-winrad);
                if( k < w )
                    continue; % Matrix is symmetric, so only compute the upper half.
                end
                
                %win_k_inds = reshape( inds( a-winrad:a+winrad, b-winrad:b+winrad ), winsize, 1 );
                win_k = reshape( I( a-winrad:a+winrad, b-winrad:b+winrad, : ), winsize, chans );
                win_kGrad = reshape( xGrad( a-winrad:a+winrad, b-winrad:b+winrad, : ), winsize, 1 );
                %win_k = Igray(win_k_inds);
                %wins(:,(c-d)+1) = win_k;

                %similarity = exp( -sum(sum((win_w - win_k).^2,2)) ./ simvar );
                if( doGrad )
                    similarity = exp( -sum(sum((win_w - win_k).^2,2).*h) ./ simvar...
                        - sum(sum(win_wGrad - win_kGrad,2).*h) ./gradVar... % Weights the central pixels higher. Adds xGradient info.
                        - ((i-a)^2 + (j-b)^2) ./ distanceVar ); % Weights distant pixels less.
                else
                    similarity = exp( -sum(sum((win_w - win_k).^2,2).*h) ./ simvar ); % Weights the central pixels higher.
                end
                %similarity = exp( -0.5*(sum(sum(abs(win_w -
                %win_k),2))/sigma - blah)^2 ); % From Bayesian Non-local means filter...doesn't work at all.
                %similarity = exp( -sum( sum((win_w - win_k).^2,2) ) ./ simvar -((i-a)^2 + (j-b)^2)/(2*(winrad)^2)); % Includes spatial term.
                similarities(pp) = similarity;
                simNdx(pp) = k;
                pp = pp + 1;
                %s( 1, c ) = similarity;
                %rowinds( 1, c ) = w;
                %colinds( 1, c ) = k;

                %c = c + 1;
            end
        end

        % Only keep top entries.
        if( pp > maxNeighbors )
            [topSims,topNdx] = weakSort(similarities, maxNeighbors, 'max');
            s(1,c:(c+maxNeighbors-1)) = topSims;
            rowinds(1,c:(c+maxNeighbors-1)) = w*ones(1,maxNeighbors);
            colinds(1,c:(c+maxNeighbors-1)) = simNdx(topNdx);
            c = c + maxNeighbors;
        else
            s(1,c:(c+pp-2)) = similarities(1:(pp-1));
            rowinds(1,c:(c+pp-2)) = w*ones(1,(pp-1));
            colinds(1,c:(c+pp-2)) = simNdx(1:(pp-1));
            c = c + pp - 1;
        end
        
        %similarity = exp( -sum( (win_w(:, ones(1,c-d)) - wins).^2 ) / simvar );
        %s( 1, d:(c-1) ) = similarity;

    end

    fprintf(1, '(%d)\n', i);
end

% Trim off rowinds, s, and colinds
lastind = size(s,2);
while( s(lastind) < 0 )
    lastind = lastind - 1;
end

if( doPatches )
    patches = patches(:,1:(c-1));
    patchNdx = patchNdx(:,1:(c-1));
else
    patches = [];
    patchNdx = [];
end

%A = sparse(rowinds(1:lastind), colinds(1:lastind), s(1:lastind), imlen, imlen);
A = sparse(rowinds(1:lastind), colinds(1:lastind), s(1:lastind), centlen, centlen);
% Only computed the upper half, so reconstruct the lower half now. Also,
% second term is (A-I) to avoid double adding the main diag.
A = A' + (A-speye(size(A)));
clear s rowinds colinds;