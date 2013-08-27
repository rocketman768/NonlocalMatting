function [ret, retNdx] = weakSort( A, k, direction )
%WEAKSORT Returns k largest or smallest elements of A.
%   A - array to be sorted
%   k - number of elements to return
%   direction - 'max' or 'min'

ret = zeros(1,k);
retNdx = zeros(1,k);
Alen = numel(A);

if( direction == 'min' )
        for i=1:k
            [minVal, minNdx] = min(A);
            ret(i) = minVal;
            retNdx(i) = minNdx;
            A(minNdx) = Inf;
        end

%     i = 0;
%     while( i < k )
%         i = i + 1;
%         minVal = A(i);
%         retNdx(i) = i;
%         for( j = i+1:Alen )
%             if( A(j) < minVal )
%                 retNdx(i) = j;
%                 i = i + 1;
%                 minVal = A(j);
%                 tmp = A(i);
%                 A(i) = A(j);
%                 A(j) = tmp;
%             end
%         end
%     end
% 
%     ret = A(1+i-k:i);
%     retNdx = retNdx(1+i-k:i);
end

if( direction == 'max' )
        for i=1:k
            [maxVal, maxNdx] = max(A);
            ret(i) = maxVal;
            retNdx(i) = maxNdx;
            A(maxNdx) = -Inf;
        end

%     i = 0;
%     while( i < k )
%         i = i + 1;
%         maxVal = A(i);
%         for( j = i+1:Alen )
%             if( A(j) > maxVal )
%                 retNdx(i) = j;
%                 i = i + 1;
%                 maxVal = A(j);
%                 tmp = A(i);
%                 A(i) = A(j);
%                 A(j) = tmp;
%             end
%         end
%     end
% 
%     ret = A(1+i-k:i);
%     retNdx = retNdx(i-k:i-1);
end

end
