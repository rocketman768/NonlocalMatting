function [ alpha, mse, W ] = matte( image, scribs, localWinRad, NLSearchRad, NLVariance, hoodSize, gt, Win, h2in )
%MATTE Gets the nonlocal matte.
%   image and scribs are required, the rest are optional.
%   localWinRad - How big the windows are when comparing image patches.
%   NLSearchRad - Radius of NL search (to limit computation).
%   NLVariance - NL Variance (bigger = more similar neighbors).
%   hoodSize - How many NL neighbors to include in linear coeff. estimation.
%   gt - The ground truth alpha matte for comparison
%   Win - The optional weight matrix so we won't have to recalculate it.
%   h2in - The h^2 parameter for Ain

if( ischar(image) )
    I = double(imread(image))./255;
else
    I = image;
end

%============YOOOHOOOOO!!==============
% Convert to grayscale
I = mean(I,3);

% Some constants.
if( ~exist('localWinRad') || isempty(localWinRad) )
    localWinRad = 3; % How big the windows are when comparing image patches.
end
if( ~exist('NLSearchRad') || isempty(NLSearchRad) )
    NLSearchRad = 10; % Radius of NL search (to limit computation).
end
if( ~exist('NLVariance') || isempty(NLVariance) )
    NLVariance = 1e0; % NL Variance (bigger = more similar neighbors).
end
if( ~exist('hoodSize') || isempty(hoodSize) )
    hoodSize = 15; % How many NL neighbors to include in linear coeff. estimation.
end
lambda = 1e-5; % Regularization param for linear coeff solution. (small)
%c = 1e5; % Regularization param for constrained alpha vals. (large)

if( ~exist('Win') || isempty(Win) || ~exist('h2in') || isempty(h2in) )
    %W = NLAdjacency(image, localWinRad, NLSearchRad, NLVariance, hoodSize);
    W = NLAdjacency(I, localWinRad, NLSearchRad, NLVariance, hoodSize);
    A = NLWeights( W );
else
    A = NLWeights( Win.^(h2in/NLVariance) );
end

if( ischar(scribs) )
    scribs = double(imread(scribs)) ./255;
end
scribs = mean(scribs,3);

[rows, cols, chans] = size(I);
imsize = rows*cols;
centrows = rows - 2*localWinRad;
centcols = cols - 2*localWinRad;
centsize = centrows * centcols;

[constr, cvals] = scribData( scribs, localWinRad );
constr = reshape( constr, centsize, 1 );
cvals = reshape( cvals, centsize, 1 );

inds = reshape([1:imsize], rows, cols);
%centinds = reshape( (1:centsize), centrows, centcols );
centToImageInds = inds( (1+localWinRad):(rows-localWinRad), (1+localWinRad):(cols-localWinRad), :);
centToImageInds = reshape( centToImageInds, centsize, 1 );
I = reshape( I, imsize, chans );

% hood = zeros( imsize, hoodSize );
% f = ones( imsize, hoodSize);

hood = zeros( centsize, hoodSize );
f = ones( centsize, hoodSize);

%for i=1:imsize
j = 1;
for i=1:centsize
    %[y, ndx] = sort( full(A(i,:)), 2, 'descend');
    %[y, hood(i,:)] = weakSort( A(:,i), hoodSize, 'max');
    [y, hood(i,:)] = weakSort( A(i,:)', hoodSize, 'max');
    %[y, hood(i,:)] = mexWeakSort( A(i,:)', hoodSize );
    %hood( i, : ) = ndx([1:hoodSize]);
    %Q = diag(y([1:hoodSize]));
    %y = y ./ sum(y(:));
    Q = diag(y);
    
    %X = I( ndx([1:hoodSize]), : );
    %X = I( hood(i,:), : );
    X = I( centToImageInds(hood(i,:)), : );
    X(:, end+1) = ones(hoodSize,1);
    
    %xprime = I( i, : );
    xprime = I( centToImageInds(i), : );
    xprime( :, end+1) = 1;
    
    % Just testing the approximate rank of X.
    % Remove these lines in practice.
    %[U,S,V] = svd(X,0);
    %rank(i) = sum( diag(S) > 0.004 ) / 4;
    
    f(i,:) = xprime * ((X'*Q*X + lambda*diag([ones(1,chans), 0])) \ (X'*Q));
    %f(i,:) = xprime * ((X'*Q*X + lambda*diag([ones(1,chans+1)])) \ (X'*Q));
    
    if( j > centsize / 100 )
        fprintf(1, '%.0f percent done\n', i/centsize*100);
        j = 1;
    else
        j = j + 1;
    end
end
clear inds I y ndx Q X xprime

fprintf(1, 'Done processing. Making matrices...\n');
%memory

%F = sparse( repmat([1:imsize]', 1, hoodSize), hood, f, imsize, imsize );
F = sparse( repmat([1:centsize]', 1, hoodSize), hood, f, centsize, centsize );
clear hood f

%% Solution Method 1. Stupid.
% C = spdiags( c * constr, 0, imsize, imsize );
% L2 = speye(imsize) - F;
% clear F;
% 
% fprintf(1, 'Getting solution...\n');
% 
% alpha = ((L2'*L2 + C) \ C) * cvals;
% alpha = reshape( alpha, rows, cols );
%% End Solution Method 1.

%% Solution Method 2. Better.
%L = speye(imsize,imsize) - F;
L = speye(centsize,centsize) - F;
clear F
L = L'*L;
% Get permutation matrix.
fprintf(1,'Getting permutation matrix...');
%P = sparse(imsize,imsize);
P = sparse(centsize,centsize);
k = sort(find(constr > 0));
j = size(k,1) + 1;
for i=1:size(k,1)
    P(i,k(i)) = 1;
end
k = sort(find(constr == 0));
for i=1:size(k,1)
    P(j-1 + i, k(i)) = 1;
end
fprintf(1,'done\n');

fprintf(1, 'Permuting laplacian...');
G = P*L*P'; % Modified laplacian.
clear L
fprintf(1, 'done\nPermuting constraints...');
vperm = P*cvals;
fprintf(1, 'done\n');

%G4 = G(j:imsize,j:imsize);
%G3 = G(j:imsize,1:j-1);
G4 = G(j:centsize,j:centsize);
G3 = G(j:centsize,1:j-1);
clear G;

ya = vperm(1:j-1);
clear vperm;

% For debugging...
%save debug.mat F G4 G3 ya;

fprintf(1, 'Solving...');
spparms('spumoni',1); % Verbose info about matrix solving.

% This works, but G4 is not usually invertible since it's basically a graph
% Laplacian, so we are really just
% asking Matlab for argmin( || G4 yb + G3 ya ||^2 ).
%yb = G4 \ (-G3*ya);

% Do LUPQ factorization to find solution since we expect that G4 can be
% decomposed into pretty sparse L and U under a permutation (P,Q).
% tic;
% [Lo,Up,p,q] = lu(G4);
% yb = q*(Up \ (Lo \ ( p*(-G3*ya) )));
% toc;
% clear Up Lo p q;

% Try Cholesky decomposition.
% tic;
% R = chol(G4);
% yb = R\(R'\(-G3*ya));
% toc;
% clear R;

% Cholesky decomposition with a sparsity inducing permutation matrix S.
% tic;
% [R,p,S] = chol(G4);
% yb = S*(R\(R'\(S'*(-G3*ya))));
% toc;
% clear R p S;

% Lower cholesky decomposition with permutation matrix S. Supposedly takes
% less memory and time. S'*G4*S = R*R'.
tic;
[R,p,S] = chol(G4, 'lower');
yb = S*(R' \ (R \ (S'*(-G3*ya))));
toc;
clear R p S;

fprintf(1, 'done\n');
alpha = P' * [ya;yb];
alpha = reshape(alpha,centrows,centcols);
%% End Solution method 2.

% Clip the matte.
alpha = min( 1, alpha );
alpha = max( 0, alpha );

%rank = reshape( rank, centrows, centcols );

% If we have a ground truth, get ssd.
if( ~exist('gt') || isempty(gt) )
    return;
end

if( ischar(gt) )
    G = double(imread(gt))./255;
else
    G = gt;
end
G = mean(G,3);

%ssd = sum(sum((alpha(localWinRad+1:rows-localWinRad,localWinRad+1:cols-localWinRad) - G(localWinRad+1:rows-localWinRad,localWinRad+1:cols-localWinRad)).^2));
mse = sum(sum((alpha - G(localWinRad+1:rows-localWinRad,localWinRad+1:cols-localWinRad)).^2)) / (centrows*centcols);
fprintf(1, 'MSE: %f\nSSD: %f\n', mse, mse*centrows*centcols);