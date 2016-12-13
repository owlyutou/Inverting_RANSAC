function [m,M] = bounds2(A,silent)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin<2)
    silent = 1;
end
[m,M] = bounds(A(:),silent);