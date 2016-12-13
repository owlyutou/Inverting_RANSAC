function [ImGridAxisPoints, GridMapMatrix] = createFirstGrid(GridStepSize,sizI2,maxScale,limitToImage)
% the function creates the first grid to be used.

maxDim2 = max(sizI2);

if ~limitToImage
    padding = 0.4*maxDim2*sqrt(maxScale); % for the case of shrinking only - at least 3/5 of transformed I1 should intersect I2    
    % padding = maxDim2*maxScale-(minOverlap*minDim2);% for the general case
else
    padding = 0;
end

ImGridAxisPoints.Xpositions = 1-padding:GridStepSize:sizI2(2)+padding;
ImGridAxisPoints.Ypositions = 1-padding:GridStepSize:sizI2(1)+padding;

%% build the grid matrix
GridMapMatrix = ones(size(ImGridAxisPoints.Ypositions,2),size(ImGridAxisPoints.Xpositions,2),'uint8');


