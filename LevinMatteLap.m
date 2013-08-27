function [L,D] = LevinMatteLap(I,epsilon,winRad) %NLMatteLap(I,A,epsilon,winRad,hoodSize)
% This is the Laplacian constructed in Levin's method, without the edge
% pixels.

if (~exist('epsilon','var'))
    epsilon=0.0000001;
end
if (isempty(epsilon))
    epsilon=0.0000001;
end
if (~exist('winRad','var'))
    winRad=1;
end
if (isempty(winRad))
    winRad=1;
end

winSize=(winRad*2+1)^2;
[rows,cols,chans]=size(I);
n=rows; m=cols;
imSize=n*m;
centrows = rows - 2*winRad;
centcols = cols - 2*winRad;
centsize = centrows * centcols;

% What are they doing here?
%consts=imerode(consts,ones(winRad*2+1));

indsM=reshape([1:imSize],rows,cols);
centInds = reshape([1:centsize],centrows,centcols);
centToImageInds = indsM( (1+winRad):(rows-winRad), (1+winRad):(cols-winRad), :);
%centToImageInds = reshape( centToImageInds, centsize, 1 );
I = reshape( I, imSize, chans );
%centI = I( centToImageInds );

% How many unlabeled pixels * winSize^2.
%tlen=sum(sum(1-consts(winRad+1:end-winRad,winRad+1:end-winRad)))*(winSize^2);
tlen = centsize*winSize^2;
%tlen = centsize*hoodSize^2;

row_inds=ones(tlen ,1);
col_inds=ones(tlen,1);
vals=zeros(tlen,1);
len=0;
for j=1:centcols
    for i=1:centrows
%for i=1:centsize
        %if (consts(i,j))
        %  continue
        %end

        %-------Levin's implementation with my notation----------
        win_inds=centInds( max(i-winRad,1):min(i+winRad,centrows),max(j-winRad,1):min(j+winRad,centcols));
        win_inds=win_inds(:);
        
        winI=I( centToImageInds(win_inds),:);
        winSize = size(winI,1);
        %winI=reshape(winI,winSize,chans);
        
        winI(:,end+1) = 1;
        %Q = diag( repmat(1/winSize, 1, winSize) );
        Q = diag( ones(1,winSize) .* 1/winSize );
        tvals = winI* ((winI'*Q*winI + epsilon/winSize*diag([ones(1,chans),0])) \ (winI'*Q));
        
        row_inds(1+len:winSize^2+len)=reshape(repmat(win_inds,1,winSize),...
             winSize^2,1);
        col_inds(1+len:winSize^2+len)=reshape(repmat(win_inds',winSize,1),...
             winSize^2,1);
        vals(1+len:winSize^2+len)=tvals(:);
        len = len + winSize^2;
        %-------------------------------------------------------
        %win_mu=mean(winI,1)';
        %win_var=inv(winI'*winI/winSize-win_mu*win_mu' +epsilon/winSize*eye(chans));

        %%win_mu = centI'*A(:,i); % Nonlocal mean. win_mu = chans x 1.
        %%win_mu = mean(winI,1)'; % Regular mean.
        %win_mu = y*winI; % Nonlocal mean.
        %winI=winI-repmat(win_mu',hoodSize,1); % Subtract the mean to get zero-mean vars.
        %%win_var = inv((winI'*winI + epsilon*eye(chans))./hoodSize); % Get inverse covariance matrix.
        %win_var = inv((winI'*Q*winI) + epsilon*eye(chans)); % Weighted inverse covariance matrix.
        
        %tvals=(1+winI*win_var*winI');

%         [y, win_inds] = weakSort( A(i,:)', hoodSize, 'max');
%         win_inds = reshape(win_inds,hoodSize,1);
%         Q = diag(y);
%         
%         winI = I( centToImageInds(win_inds), : );
%         winI(:,end+1) = 1; % Augment last column to 1.
%         tvals = winI* ((winI'*Q*winI + epsilon/hoodSize*diag([ones(1,chans),0])) \ (winI'*Q));
%         
%         row_inds(1+len:hoodSize^2+len)=reshape(repmat(win_inds,1,hoodSize),...
%             hoodSize^2,1);
%         col_inds(1+len:hoodSize^2+len)=reshape(repmat(win_inds',hoodSize,1),...
%             hoodSize^2,1);
%         vals(1+len:hoodSize^2+len)=tvals(:);
%         len=len+hoodSize^2;
    end
end

vals=vals(1:len);
row_inds=row_inds(1:len);
col_inds=col_inds(1:len);
L=sparse(row_inds,col_inds,vals,centsize,centsize);
%L=sparse(row_inds,col_inds,vals,imSize,imSize);

sumL=sum(L,2);
L=spdiags(sumL(:),0,centsize,centsize)-L;
%L=spdiags(sumL(:),0,imSize,imSize)-L;

% Normalize the matrix?
D = spdiags( spdiags(L,0), 0, centsize, centsize );
L = D \ L;
D = spdiags(D,0);

end