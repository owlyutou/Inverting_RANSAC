function res = isClockWise (x1, x2, x3)
% checks whether the points x1[1x2], x2[1x2] and x3[1x2] comprise 2 lines
% (x1->x2, x2->x3) which are clockwise relatively to each other.
% we do the comparison by checking the slopes of the 2 lines and if
% slope1 > slope2 then clockwise, otherwise - not.
% a case where slopes are the same will also return "false"
%
% assumption: input points are arrays of 1x2, where xi[1] is the x 
% coordinate, xi[2] is the y coordinate

res = (x3(2)-x2(2))*(x2(1)-x1(1)) < (x2(2)-x1(2))*(x3(1)-x2(1));
%     if (x3(2)-x2(2))*(x2(1)-x1(1)) < (x2(2)-x1(2))*(x3(1)-x2(1))
%         res = true;
%     else
%         res = false;
%     end
end