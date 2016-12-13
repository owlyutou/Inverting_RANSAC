function [drillBestXs,drillBestYs,drillBestDist,drill_best_distances,drillBestHomog] = ...
    ExpandLocally(initialExpansionRad,templateCorners,x1,x2,ImGridAxisPoints,...
    NetOfQuadrilaterals,NetOfQuadrilaterals_indices,bestInd,bestDist,GridStepSize,nInliers,byPercentileAverage)

bestImGridAxisPoints = ImGridAxisPoints;
bestNetOfQuadrilaterals = NetOfQuadrilaterals;
bestNetOfQuadrilaterals_indices = NetOfQuadrilaterals_indices;
bestbestInd = bestInd;
bestGridStepSize = GridStepSize;
drillBestDist = bestDist;
for refineRound = 1 : 400 % there is a breaking condition inside
    % % % % %
    bestConfig_indices = bestNetOfQuadrilaterals_indices(:,bestbestInd);
    bestConfig = bestNetOfQuadrilaterals(:,bestbestInd);
    
    % expand independently for factor 2 using unique
    factor = 2;
    neigh = 4; % 8
    expansionRad = 3;
    if refineRound==1
        expansionRad = initialExpansionRad;
    end
    if neigh == 8
        expansionRate = expansionRad^8;
    else % neigh == 4
        expansionRate = 5^4;
    end        
    reppedIndicesMat = kron(factor*bestConfig_indices,ones(1,expansionRate,'int16'));
    perturbMat = GetPerturbMat('int16',expansionRad); % 8 x 3^8
    reppedPerturbMat = repmat(perturbMat,[1,length(bestbestInd)]);
    
    bestNetOfQuadrilaterals_indices = reppedIndicesMat + reppedPerturbMat;
    
    vecXs = [(bestImGridAxisPoints.Xpositions - 0.5*bestGridStepSize) bestImGridAxisPoints.Xpositions(end) + 0.5*bestGridStepSize];
    xpos = zeros(1,1+2*length(bestImGridAxisPoints.Xpositions));
    xpos(1:2:end) = vecXs;
    xpos(2:2:end) = bestImGridAxisPoints.Xpositions;
    bestImGridAxisPoints.Xpositions = xpos;
    
    vecYs = [(bestImGridAxisPoints.Ypositions - 0.5*bestGridStepSize) bestImGridAxisPoints.Ypositions(end) + 0.5*bestGridStepSize];
    ypos = zeros(1,1+2*length(bestImGridAxisPoints.Ypositions));
    ypos(1:2:end) = vecYs;
    ypos(2:2:end) = bestImGridAxisPoints.Ypositions;
    bestImGridAxisPoints.Ypositions = ypos;
    
    bestNetOfQuadrilaterals = zeros(size(bestNetOfQuadrilaterals_indices),'int16');
    bestNetOfQuadrilaterals(1:2:end,:) = bestImGridAxisPoints.Xpositions(bestNetOfQuadrilaterals_indices(1:2:end,:));
    bestNetOfQuadrilaterals(2:2:end,:) = bestImGridAxisPoints.Ypositions(bestNetOfQuadrilaterals_indices(2:2:end,:));
    
    % homographies
    bestHs = SquaresToHomography_mex(templateCorners.x', templateCorners.y',double(bestNetOfQuadrilaterals(1:2:end,:)),double(bestNetOfQuadrilaterals(2:2:end,:)));
    
    % Apply homographies: Evaluate
    drill_best_distances = ApplyHomographiesOnLocations_mex(bestHs,x1(1,:),x1(2,:),x2(1,:),x2(2,:)); % mex
    drill_best_distances = sort(drill_best_distances);
    if byPercentileAverage
        drill_best_distances = bsxfun(@rdivide,cumsum(drill_best_distances,1),(1:length(x1))');
    end
    best_percentileDists = double(drill_best_distances(nInliers,:));
    
    prevBestDist = drillBestDist;
    [drillBestDist,bestbestInd] = min(best_percentileDists);
    
    newBestConfig = bestNetOfQuadrilaterals(:,bestbestInd);
    drillBestHomog = bestHs(:,bestbestInd);
    %
    % visualize the filtered grid-map and the new grid-map over the new grid:
    if 0
        figure; hold on
        title(sprintf('refinement round %d. Best score reduced from %.1f to %.1f',refineRound,prevBestDist,drillBestDist));
        axis([bx,Bx,by,By]);
        [gridXs_new,gridYs_new] = meshgrid(bestImGridAxisPoints.Xpositions,bestImGridAxisPoints.Ypositions);
        % plotting the new grid
        plot(gridXs_new,gridYs_new,'.k'); hold on;
        % plotting the old grid
        plot(gridXs_new(2:2:end,2:2:end),gridYs_new(2:2:end,2:2:end),'.k','markersize',16);
        % setting axis according to dest. image size
        for i = 1 : 4
            hh1 = plot(GT_targetCornerXs(i),GT_targetCornerYs(i),'xr','markersize',10,'linewidth',3);
            unique_corner_coords = unique(bestNetOfQuadrilaterals(2*i-1:2*i,:)','rows')';
            hh2 = plot(unique_corner_coords(1,:),unique_corner_coords(2,:),'Xb','markersize',10);
            hh3 = plot(bestConfig(2*i-1,:),bestConfig(2*i,:),'og','markersize',10);
            hh4 = plot(newBestConfig(2*i-1,:),newBestConfig(2*i,:),'or','markersize',10);
            legend([hh1,hh2,hh3,hh4],{'GT','candidates','orig','best'});
        end
    end
    % reduce local step size
    bestGridStepSize = 0.5*bestGridStepSize;
    if drillBestDist > 0.995*prevBestDist
        drillBestXs = newBestConfig(1:2:end);
        drillBestYs = newBestConfig(2:2:end);
        break
    end
    
end
