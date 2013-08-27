function [ Q ] = Q(rows, cols)
%Q Summary of this function goes here
%   Detailed explanation goes here

H = fspecial('laplacian');
[Hm,Hn] = size(H);
Hrad = 1;
N = rows*cols;

rowinds = ones(N, Hm*Hn);
colinds = ones(N,Hm*Hn);
data = zeros(N,Hm*Hn);

inds = reshape(1:N,rows,cols);

k = 1;
for i=1:rows
    for j=1:cols
        win = inds(max(1,i-Hrad):min(rows,i+Hrad),max(1,j-Hrad):min(cols,j+Hrad));
        topcutoff = max(0,1-i+Hrad);
        leftcutoff = max(0,1-j+Hrad);
        bottomcutoff = max(0,i+Hrad-rows);
        rightcutoff = max(0, j+Hrad-cols);
        [r,s] = size(win);
        data(k,1:r*s) = reshape( H((1+topcutoff):(Hm-bottomcutoff),(1+leftcutoff):(Hn-rightcutoff)), 1, r*s );
        win = reshape(win,1,r*s);
        colinds(k,1:r*s) = win;
        rowinds(k,1:r*s) = repmat([k],1,r*s);
        
        k = k + 1;
    end
end

Q = sparse(rowinds,colinds,data, N, N);

end
