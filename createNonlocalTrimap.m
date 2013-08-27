function [trimap, input, clusters] = createNonlocalTrimap(I, localWinRad, NLSearchRad, NLVariance, distanceVar, maxNeighbors)
%CREATENONLOCALTRIMAP Interactive tool to create the nonlocal trimap for matting
%  I - input normalized [0,1] image

numEigs = 10;
numClusters = 15;
[rows,cols,chans] = size(I);
centrows = rows-2*localWinRad;
centcols = cols-2*localWinRad;

%% Parse inputs
if( ~exist('maxNeighbors') || isempty(maxNeighbors) )
    maxNeighbors = (2*NLSearchRad+1)^2;
end

%% Get the nonlocal adjacency, degree, and Laplacian matrices.
A = NLAdjacency(I, localWinRad, NLSearchRad, NLVariance, maxNeighbors, 0, 0, 0, distanceVar);
D = diag(sum(A,2));
L = (D - A)' * (D - A); % Normally, L=D-A, but this seems to work well.

%% Create DD = D^(-1/2)
DD = D;
for i=1:size(D,1)
    DD(i,i) = 1/sqrt(D(i,i));
end

%% Normalize the Laplacian
M = DD*L*DD;

%% Extract smallest eigenvectors
[V,d] = eigs(M,numEigs,'SM');

%% Normalize each row of eigenvector matrix
U = V;
for i=1:size(U,1)
    U(i,:) = U(i,:) ./ norm(U(i,:));
end

%% Perform k-means clustering
[IDX,C] = kmeans(U,numClusters, 'emptyaction', 'drop');
clusters = reshape(IDX,[rows-2*localWinRad,cols-2*localWinRad]);

trimap = 0.5*ones( centrows, centcols );
input = trimap;
F = zeros(size(trimap));
B = zeros(size(trimap));

figure; imshow(I);
figure; imshow(clusters,[]);

figure;
%Foreground
[X,Y,P] = impixel( (clusters-1)./numClusters );

for i=1:length(X)
    clusterNumber = clusters(Y(i),X(i));
    F( find(clusters==clusterNumber) ) = 1;
    input(Y(i),X(i)) = 1;
end

%Background
[X,Y,P] = impixel( (clusters-1)./numClusters );

for i=1:length(X)
    clusterNumber = clusters(Y(i),X(i));
    B( find(clusters==clusterNumber) ) = 1;
    input(Y(i),X(i)) = 0;
end



%% Erode the edges.
SE = ones(2*localWinRad+1);
F = imerode(F, SE);
B = imerode(B, SE);

%% Create Trimap.
trimap( find(F == 1) ) = 1;
trimap( find(B == 1) ) = 0;

imshow(trimap);

clusters = (clusters-1)./numClusters;

end