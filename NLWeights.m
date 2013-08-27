function [L] = NLWeights( A )
% This is not really a laplacian since
% that would be I - A instead of I + A.
% This simply inserts a 1 on the diagonal, meaning that pixel i is exactly
% similar to pixel i, then row normalizes so that the rows of L sum to 1.
% In other words, this is a normalized kernel matrix.

%L = speye(size(A)) + A;
[rows,cols] = size(A);

% This is the random walk Laplacian.
L = spdiags( sum(A,2), 0, rows, cols) \ A;

% This is a normalized Laplacian.
%D_1_2 = spdiags( sqrt(sum(L,2)), 0, rows, cols);
%L = D_1_2 \ L / D_1_2;