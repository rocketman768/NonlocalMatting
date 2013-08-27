function [L] = NLLaplacian( A )
% Takes an adacency matrix A and outputs the Laplacian L.

% NOTE TO SELF: since the denoising procedure also utilizes the underlying
% pixel, perhaps I should put A = A + I here.
A = A + speye(size(A));
[rows,cols] = size(A);
D = spdiags( sum(A,2), 0, rows, cols);
L = D - A;

% This is the random walk Laplacian: L = I - D^{-1}A
%L = D \ L;

% This is a normalized Laplacian.
%D_1_2 = spdiags( sqrt(sum(L,2)), 0, rows, cols);
%L = D_1_2 \ L / D_1_2;