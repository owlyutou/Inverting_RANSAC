function [NetOfQuadrilaterals,NetOfQuadrilaterals_indeces,CurrentNetIndex,OutOfMem] = ...
    build_net_shrinking_homographies(sizI1,sizI2,scaleFact,ImGrid,CornerValidMaps,GridStepSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About this function:
% This function builds an alfa-cover Net of homogeneous transforms in the
% destination image.
% For the Net to be an alfa-cover, we'll divide the destination image into
% a grid with a block size of floor([alfa x alfa]). This
% way, each transformation which maps the original template to the
% destination image will have a transformation that is up to alfa far away.
%
% The transformations are viewed in the following way:
%   1. The template has its 4 corner coordinates.
%   2. The transformed template in the destination image plane has new 4
%      corners.
%   3. The form of the resulting shape after the transformation is a
%      quadrilateral.
%   4. At first implementation phase, we'll collect all the quadrilaterals
%      in the destination image.
%        4.1 There will be some restrictions on valid quadrilaterals (not
%            to get extreme shapes).
%        4.2 The restrictions will be parametrical.
%        4.3 The main restriction: the 4 corner points of the quadrilateral
%            must lay on the grid points of the destination image.
%   5. At second implementation phase, we'll find the Transformation
%      Matrices that bring the template rectangle to each of the
%      quadrilaterals found before. Eventually having one transformation
%      matrix per valid quadrilateral.
%
% The scheme of constructing the quadrilaterals:
%   1. Set the first point of the quadrilateral (can be any point of the
%      grid).
%   2. Set the second point of the quadrilateral:
%        2.1 Must be on a grid point of the destination image.
%        2.2 Must be no more than "L_X2_X1_MAX_DIST" (tunable parameter)
%            pixels far from the first point.
%   3. Set the third point of the quadrilateral:
%        3.1 Must be on a grid point of the destination image.
%        3.2 Must be no more than "L_X3_X2_MAX_DIST" (tunable parameter)
%            pixels far from the second point.
%        3.3 Must be no more than "L_X3_X1_MAX_DIST" (tunable parameter)
%            pixels far from the first point.
%        3.4 Must be no LESS than "L_X3_X1_MIN_DIST" (tunable parameter)
%            pixels far from the first point.
%   4. Set the 4th point of the quadrilateral:
%        4.1 Must be on a grid point of the destination image.
%        4.2 Must be no more than "L_X4_X1_MAX_DIST" (tunable parameter)
%            pixels far from the first point.
%        4.3 Must be no more than "L_X4_X3_MAX_DIST" (tunable parameter)
%            pixels far from the third point.
%        4.4 Must be no more than "L_X4_X2_MAX_DIST" (tunable parameter)
%            pixels far from the second point.
%        4.5 Must be no LESS than "L_X4_X2_MIN_DIST" (tunable parameter)
%            pixels far from the second point.
% In overall, we have 8 degrees of freedom for which we got 8 parameters.
%
% The inputs are:
%   i_N1 - this is the template size in pixels. assuming the template is
%          of size [i_N1 x i_N1] pixels.
%   ImGridAxisPoints.Xpositions - vector of 1xn. x values of the grid
%                                 points.
%   ImGridAxisPoints.Ypositions - vector of 1xm. y values of the grid
%                                 points.
%   CornerValidMaps(i).GridMapMatrix - matrix of size mxn. in each cell of the matrix there is
%                   a "1" for a valid grid point. since this is the first
%                   grid generated - all points are valid. thus - this is a
%                   ones matrix of size mxn. type: uint8
%   GridStepSize - vertical and horizontal distance (in destination image
%                  pixels) between 2 adjescent grid points. (size: 1x1,
%                  double)
% The outputs are:
%   NetOfQuadrilaterals - [8xm] matrix of m quadrilaterals. each
%                         quadrilateral is defined by 4 2D points (hense
%                         the 8). The 2D points are in terms of the
%                         destination image plane (1x1..dest_im_size).
%   NetOfQuadrilaterals_indeces - same as NetOfQuadrilaterals in size and
%                                 meaning, but the 2D points are presented
%                                 as the index numbers of ImGridAxisPoints:
%                                 the x points are indeces of the vector
%                                 ImGridAxisPoints.Xpositions
%                                 the y points are indeces of the vector
%                                 ImGridAxisPoints.Ypositions
%   OutOfMem - indicates if the preallocated space for NetOfQuadrilaterals
%              was sufficient (0) or not (1). If wasn't suffitient - the
%              returned Net is at the allocated size, but doesn't cover all
%              the possible quadrilaterals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if coder.target('MATLAB')
    % Code for MATLAB execution, ignored by Matlab coder
    if exist('build_net_shrinking_homographies_mex','file')==3
        [NetOfQuadrilaterals,NetOfQuadrilaterals_indeces,CurrentNetIndex,OutOfMem] = ...
            build_net_shrinking_homographies_mex(sizI1,sizI2,scaleFact,ImGrid,CornerValidMaps,GridStepSize);
        return
    else
        warning('build_net_shrinking_homographies:NoMex','!!! no mex file found, running matlab version is ~x100 times slower!!!')
    end
end

% initialize parameters:

h1 = sizI1(1);
w1 = sizI1(2);

minVertDist = h1/scaleFact;
maxVertDist = h1*scaleFact;
minHorizDist = w1/scaleFact;
maxHorizDist = w1*scaleFact;
minDiagDist = norm(sizI1)/scaleFact;
maxDiagDist = norm(sizI1)*scaleFact;
minDiagRatio = 1/scaleFact;
maxDiagRatio = scaleFact;

minRatio_vert_hor = (h1/w1)/scaleFact^1.2;
maxRatio_vert_hor = (h1/w1)*scaleFact^1.2;
minRatio_opposite = 1/scaleFact;
maxRatio_opposite = 1*scaleFact;

colinearityFact = 1.05;

% area considerations
areaI2square = prod(sizI2)^2;
areaExpansionFactSqr = 1.1^2;

% indicate if memory allocated was enough(0) or not(1)
OutOfMem = 0;

maxStepsX = ceil(w1*scaleFact/GridStepSize);
maxStepsY = ceil(h1*scaleFact/GridStepSize);

% make a non dynamic allocation (for debug):
CurrentNetSize = 10e6;
NetOfQuadrilaterals_indeces = zeros(8,CurrentNetSize,'int16');
NetOfQuadrilaterals = zeros(8,CurrentNetSize,'int16');

CurrentNetIndex = 0;

%% build the net
% go over the possible first coordinate
% possible coordinates: those that fall on a grid point in the output plane.

% -- P1 --
for idx_x1 = 1:length(ImGrid.Xpositions)
    x1Pos_x = ImGrid.Xpositions(idx_x1);
    % fprintf('x1Pos_x = %.0f  out of the range: %.0f : %.0f : %.0f\n',x1Pos_x,ImGrid.Xpositions(1),GridStepSize-1,ImGrid.Xpositions(end));
%     fprintf('x1Ind = %d  out of: %d\n',idx_x1,length(ImGrid.Xpositions));
    for idx_y1 = 1:length(ImGrid.Ypositions)
        if (~CornerValidMaps(1).GridMapMatrix(idx_y1,idx_x1)) % check if point1 is available in the grid map
            continue
        end
        x1Pos_y = ImGrid.Ypositions(idx_y1);
        % -- P2 --
        for idx_x2 = max(1,idx_x1-maxStepsX):min(length(ImGrid.Xpositions),idx_x1+maxStepsX)
            x2Pos_x = ImGrid.Xpositions(idx_x2);
            for idx_y2 = max(1,idx_y1-maxStepsY):min(length(ImGrid.Ypositions),idx_y1+maxStepsY)
                if (~CornerValidMaps(2).GridMapMatrix(idx_y2,idx_x2)) % check if point2 is available in the grid map
                    continue
                end
                x2Pos_y = ImGrid.Ypositions(idx_y2);
                ed_12 = euclidDist([x1Pos_x,x1Pos_y],[x2Pos_x,x2Pos_y]);
                if ed_12 < minVertDist || ed_12 > maxVertDist
                    continue;
                end
                % -- P3 --
                for idx_x3 = max(1,idx_x2-maxStepsX):min(length(ImGrid.Xpositions),idx_x2+maxStepsX)
                    x3Pos_x = ImGrid.Xpositions(idx_x3);
                    for idx_y3 = max(1,idx_y2-maxStepsY):min(length(ImGrid.Ypositions),idx_y2+maxStepsY)
                        if (~CornerValidMaps(3).GridMapMatrix(idx_y3,idx_x3)) % check if point3 is available in the grid map
                            continue
                        end
                        x3Pos_y = ImGrid.Ypositions(idx_y3);
                        % first check that the 3 points construct 2
                        % clockwise lines (going from x1 to x3)
                        if ~isClockWise([x1Pos_x, x1Pos_y],[x2Pos_x, x2Pos_y],[x3Pos_x, x3Pos_y])
                            continue
                        end
                        ed_13 = euclidDist([x1Pos_x,x1Pos_y],[x3Pos_x,x3Pos_y]);
                        if ed_13 < minDiagDist || ed_13 > maxDiagDist
                            continue
                        end
                        ed_23 = euclidDist([x2Pos_x,x2Pos_y],[x3Pos_x,x3Pos_y]);
                        if ed_23 < minHorizDist || ed_23 > maxHorizDist
                            continue
                        end
                        if (ed_12/ed_23) < minRatio_vert_hor || (ed_12/ed_23) > maxRatio_vert_hor
                            continue
                        end

                        % -- P4 --
                        for idx_x4 = max(1,idx_x3-maxStepsX):min(length(ImGrid.Xpositions),idx_x3+maxStepsX)
                            x4Pos_x = ImGrid.Xpositions(idx_x4);
                            for idx_y4 = max(1,idx_y3-maxStepsY):min(length(ImGrid.Ypositions),idx_y3+maxStepsY)
                                if (~CornerValidMaps(4).GridMapMatrix(idx_y4,idx_x4)) % check if point4 is available in the grid map
                                    continue
                                end
                                x4Pos_y = ImGrid.Ypositions(idx_y4);
                                % first check that the points construct a
                                % clockwise quadrilateral
                                if ~isClockWise([x2Pos_x, x2Pos_y],[x3Pos_x, x3Pos_y],[x4Pos_x, x4Pos_y]) ...
                                        || ~isClockWise([x3Pos_x, x3Pos_y],[x4Pos_x, x4Pos_y],[x1Pos_x, x1Pos_y]) ...
                                        || ~isClockWise([x4Pos_x, x4Pos_y],[x1Pos_x, x1Pos_y],[x2Pos_x, x2Pos_y])
                                    continue
                                end
                                ed_14 = euclidDist([x1Pos_x,x1Pos_y],[x4Pos_x,x4Pos_y]);
                                if ed_14 < minHorizDist || ed_14 > maxHorizDist
                                    continue
                                end
                                ed_24 = euclidDist([x2Pos_x,x2Pos_y],[x4Pos_x,x4Pos_y]);
                                if ed_24 < minDiagDist || ed_24 > maxDiagDist
                                    continue
                                end
                                ed_34 = euclidDist([x3Pos_x,x3Pos_y],[x4Pos_x,x4Pos_y]);
                                if ed_34 < minVertDist || ed_34 > maxVertDist
                                    continue;
                                end
                                diagRatio = ed_13/ed_24;
                                if diagRatio < minDiagRatio || diagRatio > maxDiagRatio
                                    continue
                                end
                                if (ed_12/ed_34) < minRatio_opposite || (ed_12/ed_34) > maxRatio_opposite
                                    continue
                                end
                                if (ed_12/ed_14) < minRatio_vert_hor || (ed_12/ed_14) > maxRatio_vert_hor
                                    continue
                                end
                                if (ed_34/ed_23) < minRatio_vert_hor || (ed_34/ed_23) > maxRatio_vert_hor
                                    continue
                                end
                                if (ed_23/ed_14) < minRatio_opposite || (ed_23/ed_14) > maxRatio_opposite
                                    continue
                                end
                                if (ed_34/ed_14) < minRatio_vert_hor || (ed_34/ed_14) > maxRatio_vert_hor
                                    continue
                                end
                                % colinearity
                                if ed_12 + ed_23 <= colinearityFact * ed_13
                                    continue
                                end
                                % colinearity
                                if ed_34 + ed_14 <= colinearityFact * ed_13
                                    continue
                                end
                                % colinearity
                                if ed_14 + ed_12 <= colinearityFact * ed_24
                                    continue
                                end
                                % colinearity
                                if ed_23 + ed_34 <= colinearityFact * ed_24
                                    continue
                                end

                                % limit area expansion, using Brahmagupta's formula (known from ~600 CE :-)
                                semiPerim = 0.5*(ed_12+ed_23+ed_34+ed_14);
                                targetAreaSquare = (semiPerim-ed_12)*(semiPerim-ed_23)*(semiPerim-ed_34)*(semiPerim-ed_14);

                                if targetAreaSquare > areaI2square*areaExpansionFactSqr
                                    continue
                                end


                                % if came to here - we have a valid quadrilateral
                                % we'll save it
                                % first: check that we have enough space:
                                CurrentNetIndex = CurrentNetIndex + 1;
                                if CurrentNetIndex > CurrentNetSize % shouldn't happen
                                    CurrentNetIndex = CurrentNetIndex - 1;
                                    % leave the current net, and report out of mem
                                    OutOfMem = 1;
                                    return;
                                end
                                NetOfQuadrilaterals(:,CurrentNetIndex) = int16(0.5 + ...
                                    [x1Pos_x; x1Pos_y; x2Pos_x; x2Pos_y;...
                                    x3Pos_x; x3Pos_y; x4Pos_x; x4Pos_y]);
                                NetOfQuadrilaterals_indeces(:,CurrentNetIndex) = int16(...
                                    [idx_x1; idx_y1; idx_x2; idx_y2;...
                                     idx_x3; idx_y3; idx_x4; idx_y4]);
                            end
                        end
                    end
                end
            end %for x2Pos_y
        end %for x2Pos_x
    end %for x1Pos_y
end %for x1Pos_x

% the cropping of the buffer will be done outsize of this function

fprintf('Net size: %.0f\n',CurrentNetIndex);


end % of function build_net





