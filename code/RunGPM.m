function GPMresult = RunGPM(GPMparams,x1,x2,OutputDir,...
        GT_targetCornerXs,GT_targetCornerYs,targetCorners,templateCorners,orig_GT_distances,I1,I2,USAC_distances)
% GPMresult contains: 
% [status,drillBestXs,drillBestYs,H_drill,distances_drill,drillBestDist,percentileToUse,nInliers,lowerBound,BnBlevels,GPM_result_fig]

%% GPM 0 - some local params

GPMresult = [];

doPlots = 0;
numMatches = size(x1,2);

doIRE = @(needToEvaluateIRE,level)((needToEvaluateIRE && (level <=6)) || (level <=4)); %(level==1);
guarantyFactor = 1/sqrt(2); % relation between grid step and error guarantee
maxBnBlevels = 20;
debugVis = false; % true; %

% image sizes
[h1,w1,d1] = size(I1);
[h2,w2,d2] = size(I2);

% unpack external params
byPercentileAverage = GPMparams.byPercentileAverage; % 1 means by-average
if isfield(GPMparams,'robustPercentile')
robustPercentile = GPMparams.robustPercentile; % value in [0,1]
else
robustPercentile = 1; % value in [0,1]
end
scaleFact = GPMparams.scaleFact; % maximum scale

% just for reference (not really used)
if byPercentileAverage
    GT_distances = sort(orig_GT_distances);
    GT_distances = bsxfun(@rdivide,cumsum(GT_distances,1),(1:length(x1))');
    if exist('USAC_distances','var');
        used_USAC_distances = sort(USAC_distances);
        used_USAC_distances = bsxfun(@rdivide,cumsum(used_USAC_distances,1),(1:length(x1))');
    end
end

%% GPM 1 - check that GT scale is in bounds (including a figure)


% measure overlap
[xs_inter, ys_inter] = polybool('intersection',GT_targetCornerXs',GT_targetCornerYs',targetCorners.x,targetCorners.y);
GT_overlap = polyarea(xs_inter,ys_inter)/(h2*w2);
% measure scale
[indsSource,indsTarget] = meshgrid(1:4,1:4);
sourceDists = sqrt( (templateCorners.x(indsSource) - templateCorners.x(indsTarget)).^2 + ...
    (templateCorners.y(indsSource) - templateCorners.y(indsTarget)).^2   );
gtDists =     sqrt( (GT_targetCornerXs(indsSource) - GT_targetCornerXs(indsTarget)).^2 + ...
    (GT_targetCornerYs(indsSource) - GT_targetCornerYs(indsTarget)).^2   );
[GTminScale,GTmaxScale] = bounds2(gtDists./sourceDists);


% the Figure
if doPlots
    gtScaleFig = figure; hold on;
    plot([targetCorners.x targetCorners.x(1)],[targetCorners.y targetCorners.y(1)],'m','linewidth',5);
    plot([GT_targetCornerXs' GT_targetCornerXs(1)],[GT_targetCornerYs' GT_targetCornerYs(1)],'g','linewidth',5)
    plot([xs_inter xs_inter(1)],[ys_inter ys_inter(1)],':b','linewidth',3)
    xBounds = bounds2([targetCorners.x GT_targetCornerXs']) + [-10,10];
    yBounds = bounds2([targetCorners.y GT_targetCornerYs']) + [-10,10];
    axis([xBounds(1),xBounds(2),yBounds(1),yBounds(2)]);
    title(sprintf('GT minScale: %.2f, maxScale: %.2f, overlap: %.2f',GTminScale,GTmaxScale,GT_overlap));
    legend('I_2','GT','intersection','location','best');
    title('Checking GT scale')
    axis equal; axis tight
end

% abort this example if scale is not in bounds
if ~(max(GTmaxScale,1/GTminScale) < scaleFact)
    warning('!!! GT transformation (%.1f) is out of our scale bounds (%.1f) - SKIPPING EXAMPLE!!!\n',scaleFact,max(GTmaxScale,1/GTminScale));
    GPMresult.status = 'GT transformation is out of our scale bounds';
    return
end



%% GPM 2 - Create initial net

% 0] params

StaticAllocation = 1;
useMex = 1;
GridStepSize = min([h1,w1])/3; % was 3;
minOverlap = 0.9; % obsolete now
limitToImage = 0;

paramsStr = sprintf('_%.2f_%d_%.1f_',scaleFact,GridStepSize,minOverlap);

CornerValidMaps = [];
needToEvaluateIRE = 1;
firstSuccessfullIRElevel = 999;

% 1] the points grid

[ImGridAxisPoints, GridMapMatrix] = createFirstGrid(GridStepSize,[h2,w2],scaleFact,limitToImage);
for i = 1 : 4
    CornerValidMaps(i).GridMapMatrix = GridMapMatrix;
end

% 2] the net

[NetOfQuadrilaterals,NetOfQuadrilaterals_indices,CurrentNetIndex,OutOfMem] = ...
    build_net_shrinking_homographies([h1,w1],[h2,w2],scaleFact,ImGridAxisPoints,CornerValidMaps,GridStepSize);

assert(~OutOfMem)
while CurrentNetIndex > 0.8*10^6 % was 0.8*10^6
    fprintf('too large grid... GridStepSize: %.1f, CurrentNetIndex: %d\n',GridStepSize,CurrentNetIndex);
    GridStepSize = 1.05*GridStepSize;
    [ImGridAxisPoints, GridMapMatrix] = createFirstGrid(GridStepSize,[h2,w2],scaleFact,limitToImage);
    for i = 1 : 4
        CornerValidMaps(i).GridMapMatrix = GridMapMatrix;
    end
    %             [NetOfQuadrilaterals,NetOfQuadrilaterals_indices,CurrentNetIndex,OutOfMem] = ...
    %                 build_net_homography_mex(size(I1),scaleFact,ImGridAxisPoints,CornerValidMaps,GridStepSize);
    [NetOfQuadrilaterals,NetOfQuadrilaterals_indices,CurrentNetIndex,OutOfMem] = ...
        build_net_shrinking_homographies_mex([h1,w1],[h2,w2],scaleFact,ImGridAxisPoints,CornerValidMaps,GridStepSize);
end
% warning('Reminder (SIMON): NetOfQuadrilaterals is currently int16 instead of double to save space...');

% crop the buffer, keep only valid results:
NetOfQuadrilaterals = NetOfQuadrilaterals(:,1:CurrentNetIndex);
NetOfQuadrilaterals_indices = NetOfQuadrilaterals_indices(:,1:CurrentNetIndex);


%% 5] BRANCH and BOUND

% for visualizing only:
[bx,Bx] = bounds2(ImGridAxisPoints.Xpositions);
[by,By] = bounds2(ImGridAxisPoints.Ypositions);

% initializations
levelStats = [];
drillBestDists = [];
lowerBound = 0; % initial lower bound
percentileToUse = -99; % invalid value for initialization

% main loop
for level = 1 : maxBnBlevels
    
    %         close all
    
    goodIREevaluation = 0;
    
    
    %% BnB 1 - find a homography matrix for each quadrilateral.


    H = SquaresToHomography(templateCorners.x', templateCorners.y',double(NetOfQuadrilaterals(1:2:end,:)),double(NetOfQuadrilaterals(2:2:end,:)));
    
    
    %% BnB 2 - Evaluate the homographies to get the distances matrix
    tic
    if size(H,2) > 10^6
        distances = ApplyHomographiesOnLocations_int16_mex(H,x1(1,:),x1(2,:),x2(1,:),x2(2,:)); % mex
    else
        distances = ApplyHomographiesOnLocations_mex(H,x1(1,:),x1(2,:),x2(1,:),x2(2,:)); % mex
    end
    applyTime = toc;
    
    
    [distances,sortIdxs] = sort(distances); % sorting columns (errors per transformation)
    
    orig_distances = distances; % this is to diferentiate from the average version, that might be computed next
    if byPercentileAverage
        distances = bsxfun(@rdivide,cumsum(distances,1),cast((1:length(x1))',class(distances)));
    end
    
    
    %% BnB 3 - Drill for IRE
    
    if doIRE(needToEvaluateIRE,level)
        % 1] For each percentile from a range -> drill in to improve 'rOpt'
        best_distances_per_percentile = 99999999999*ones(numMatches,1);
        initialExpansionRad = 3;
        origStep = 4;
        if level == 1
            pctilerange = origStep:origStep:98;
        else
            step = origStep/2^(level-1);
            pctilerange = BoundBy(percentileToUse + (-10*step:step:10*step),2,98);
        end
        
        for GTpercentileForIRE = pctilerange %27%[24,30]
            %
            nInliersGT = round(GTpercentileForIRE*size(orig_distances,1)/100);
            percentileDists = double(orig_distances(nInliersGT,:));
            % bests for current percentile
            [bestDistForPercentile,bestIndforPercentile] = min(percentileDists);
            %
            
            localByPercentileAverage = 0; % we don't want the IRE to use averages
            [~,~,~,tmp_drill_best_distances] = ...
                ExpandLocally(initialExpansionRad,templateCorners,x1,x2,ImGridAxisPoints,NetOfQuadrilaterals,...
                NetOfQuadrilaterals_indices,bestIndforPercentile,bestDistForPercentile,GridStepSize,nInliersGT,localByPercentileAverage);
            best_distances_per_percentile = min(best_distances_per_percentile,double(min(tmp_drill_best_distances,[],2)));
            %             plot(double(min(tmp_drill_best_distances,[],2))','.r')
        end
        
        % 2] automatic inlier-rate computation
        bestVals = best_distances_per_percentile;
        nExpansions = sum(bsxfun(@lt,orig_distances,bestVals+guarantyFactor*GridStepSize),2);
        validPercentileInds = find((bestVals-guarantyFactor*GridStepSize) < 5);
        minExpansionsIndex = find(nExpansions(validPercentileInds)==min(nExpansions(validPercentileInds)),1,'last');
        chosenPercentile = validPercentileInds(minExpansionsIndex)/numMatches;
        
        badPercentile = 0;
        if isempty(chosenPercentile)
            chosenPercentile = 0.5;
            badPercentile = 1;
        end
        
        if needToEvaluateIRE || (level <=firstSuccessfullIRElevel+1)
            percentileToUse = robustPercentile*100*chosenPercentile;
        end
        
        if ~badPercentile && minExpansionsIndex < 0.99*(validPercentileInds(end))
            % updating p*
            if needToEvaluateIRE
                firstSuccessfullIRElevel = level;
            end
            fprintf('good IRE: found p* to be %.2f\n',100*chosenPercentile);
            needToEvaluateIRE = needToEvaluateIRE && 0;
            goodIREevaluation = 1;
        else
            fprintf('bad IRE: found p* to be %.2f\n',100*chosenPercentile);
            needToEvaluateIRE = needToEvaluateIRE && 1;
        end
    end
    
    %% BnB 4 - Some important transformations: "Best" as well as "GT", "NN"
    
    % A] NN
    
    gridCornersX = double(NetOfQuadrilaterals(1:2:7,:));
    gridCornersY = double(NetOfQuadrilaterals(2:2:8,:));
    
    distsX2 = (bsxfun(@minus,gridCornersX,GT_targetCornerXs)).^2;
    distsY2 = (bsxfun(@minus,gridCornersY,GT_targetCornerYs)).^2;
    sqrDists = sum(distsX2+distsY2);
    
    [nnDist,nnInd] = min(sqrDists);
    
    nnCornersX = gridCornersX(:,nnInd);
    nnCornersY = gridCornersY(:,nnInd);
    nnHomog = H(:,nnInd);
    
    NN_distances = ApplyHomographiesOnLocations_mex(nnHomog,x1(1,:),x1(2,:),x2(1,:),x2(2,:)); % mex
    
    if byPercentileAverage
        NN_distances = sort(NN_distances);
        NN_distances = bsxfun(@rdivide,cumsum(NN_distances,1),(1:length(x1))');
    end
    
    % sanity - check NN correctness by snaping GT to grid (and measuring true scale and overlap)
    snappedGTxs = zeros(4,1); snappedGTys = zeros(4,1);
    for i = 1 : 4
        [closestAbsDistX,closestIndX] = min(abs(ImGridAxisPoints.Xpositions-GT_targetCornerXs(i)));
        snappedGTxs(i) = ImGridAxisPoints.Xpositions(closestIndX);
        [closestAbsDistY,closestIndY] = min(abs(ImGridAxisPoints.Ypositions-GT_targetCornerYs(i)));
        snappedGTys(i) = ImGridAxisPoints.Ypositions(closestIndY);
    end
    
    % this is the case that the nearest in the net isn't the GT 'snapped' to the grid - find out why
    if ~(all(abs(snappedGTxs-nnCornersX)<=1) && all(abs(snappedGTys-nnCornersY)<=1))
        %%
        if 0
            %%
            build_net_shrinking_homographies_debug([h1,w1],[h2,w2],scaleFact,...
                ImGridAxisPoints,CornerValidMaps,GridStepSize,snappedGTxs,snappedGTys);
        end
        %%
        fprintf('nearest ~= snapped\n');
        if doPlots && debugVis
            nearest2snappedMismatchFig = figure(11); clf
            [gridXs,gridYs] = meshgrid(ImGridAxisPoints.Xpositions,ImGridAxisPoints.Ypositions);
            plot(gridXs,gridYs,'.k'); hold on;
            axis([bx,Bx,by,By]);
            
            hold on;
            ps = [];
            ps(end+1) = plot(GT_targetCornerXs([1:4 1]),GT_targetCornerYs([1:4 1]),'r','linewidth',2);
            ps(end+1) = plot(nnCornersX([1:4 1]),nnCornersY([1:4 1]),'c','linewidth',2);
            ps(end+1) = plot(targetCorners.x([1:4 1]),targetCorners.y([1:4 1]),'k:','linewidth',1);
            ps(end+1) = plot(snappedGTxs([1:4 1]),snappedGTys([1:4 1]),'y','linewidth',2);
            legend(ps,{'GT','Nearest','I2','snapped'})
            validAreaRatio = polyarea(snappedGTxs([1:4 1]),snappedGTys([1:4 1]))/(h1*w1);
            validAreaRatio = polyarea(nnCornersX([1:4 1]),nnCornersY([1:4 1]))/(h1*w1);
            if ~quiet
                keyboard
            end
            
            config = zeros(8,1);
            config(1:2:end) = snappedGTxs;
            config(2:2:end) = snappedGTys;
            sum(ismember(NetOfQuadrilaterals',config','rows'));
        end
    end
    
    % B] GT - distances are computed at beginning (since they aren't affected by BnB)

    % C] BEST (using our estimated 'percentileToUse')
    nInliers = max(1,round(percentileToUse*size(distances,1)/100));
    percentileDists = double(distances(nInliers,:));
    
    [bestDist,bestInd] = min(percentileDists);
    bestCornersX = gridCornersX(:,bestInd);
    bestCornersY = gridCornersY(:,bestInd);
    % specificErrors = distances(:,bestInd);
    
    %             NN_percentile = prctile(NN_distances,GTpercentile);
    %             GT_percentile = prctile(GT_distances,GTpercentile);
    NN_bestDist = NN_distances(nInliers);
    GT_bestDist = GT_distances(nInliers);
    
    
    %% BnB 5 - EXPAND best configuration locally to improve BnB threshold
    initialExpansionRad = 3;
    [drillBestXs,drillBestYs,drillBestDist,drill_best_distances,H_drill] = ...
        ExpandLocally(initialExpansionRad,templateCorners,x1,x2,ImGridAxisPoints,NetOfQuadrilaterals,...
        NetOfQuadrilaterals_indices,bestInd,bestDist,GridStepSize,nInliers,byPercentileAverage);

    if drillBestDist < bestDist
        fprintf(' :-) drilling helped - reduced error from %.1f to %.1f\n',bestDist,drillBestDist);
    end

        %% BnB 5.5 - Run local optimization (reweighted LS) to improve BnB threshold
        
        % optimization improvement
        H_opt = robustfit_Homog(H_drill,x1,x2,nInliers);
        OPTdistances_drill = ApplyHomographiesOnLocations_mex(H_opt,x1(1,:),x1(2,:),x2(1,:),x2(2,:)); % mex
        sOPTDistances_drill = sort(OPTdistances_drill);
        if byPercentileAverage
            sOPTDistances_drill = bsxfun(@rdivide,cumsum(sOPTDistances_drill,1),(1:length(x1))');
        end
        % switch to optimized, if better in terms of the upper bound it achieves
        if sOPTDistances_drill(nInliers) <= drillBestDist
            fprintf(' :-) :-) :-) reweighted LS helped - reduced error from %.1f to %.1f\n',drillBestDist,sOPTDistances_drill(nInliers));
            drillBestDist = sOPTDistances_drill(nInliers);
            H_drill = H_opt;
            zVec = [-1,-1,-1,-1];
            [~,drillBestXs,drillBestYs] = ApplyHomographiesOnLocationsFull_mex(...
                reshape(H_opt',9,1),templateCorners.x,templateCorners.y,zVec,zVec); % mex
        end

    %% BnB 6 - Compute lower-bound
    if goodIREevaluation
        lowerBound = max(0,bestDist - guarantyFactor*GridStepSize);
    else
        lowerBound = max(lowerBound,bestDist - guarantyFactor*GridStepSize);
    end
    
    if lowerBound > min(drillBestDist,GT_bestDist)
        warning('lowerbound doesn''t reach below the GT error (probably due to reevaluation of IRE)');
    end
    
    
    %% BnB 7 - Calculate the final BnB threshold and the resulting lower bound
    thresh = min(bestDist,drillBestDist) + guarantyFactor*GridStepSize;
    
    origThresh = thresh;
    
    assert(thresh > drillBestDist)
    
    tooManyExpansions = 1;
    wIters = 0;
    while(tooManyExpansions)
        % % debug
        wIters = wIters + 1;
        if wIters>50
            keyboard
        end
        % %
        % check the expansion size and refine search if too large
        goodInds = find(percentileDists < thresh);
        if length(goodInds)<=2000
            
            FilteredNetOfQuadrilaterals_indeces = NetOfQuadrilaterals_indices(:,goodInds);
            
            factor = 2;
            neigh = 4; %8
            expansionRad = 3;
            if neigh == 8
                expansionRate = expansionRad^8;
            else % neigh == 4
                expansionRate = 5^4;
            end
            reppedIndicesMat = kron(factor*FilteredNetOfQuadrilaterals_indeces,ones(1,expansionRate,'int16'));
            perturbMat = GetPerturbMat('int16',expansionRad,neigh); % 8 x 3^8
            reppedPerturbMat = repmat(perturbMat,[1,length(goodInds)]);
            
            expandedIndices = reppedIndicesMat + reppedPerturbMat;
            tmpNetOfQuadrilaterals_indices = unique(expandedIndices','rows')';
        end
        if length(goodInds)>2000 || size(tmpNetOfQuadrilaterals_indices,2) > 0.5*10^6
            %                 if size(tmpNetOfQuadrilaterals_indices,2) > 0.5*10^6
            tooManyExpansions = 1;
            prevThresh = thresh;
            thresh = 0.95*thresh;
            fprintf('EXPLOSION at level %d: reducing threshold from %.1f to %.1f\n',level,prevThresh,thresh);
            
            %                 % repeat drilling with larger neighborhood
            %                 initialExpansionRad = 5;
            %                 [drillBestXs,drillBestYs,drillBestDist] = ...
            %                     ExpandLocally(initialExpansionRad,templateCorners,x1,x2,ImGridAxisPoints,NetOfQuadrilaterals,...
            %                     NetOfQuadrilaterals_indices,bestInd,bestDist,GridStepSize,percentileToUse,byPercentileAverage);
            %                 prevThresh = thresh;
            %                 thresh = min(bestDist,drillBestDist) + guarantyFactor*GridStepSize;
            %                 warning('Performed wider drilling. Improved threshold from %.1f to %.1f',prevThresh,thresh);
        else
            tooManyExpansions = 0;
        end
    end
    
    % stats
    levelStats(level).thresh = thresh; %#ok<*AGROW,*SAGROW>
    levelStats(level).bestDist = bestDist;
    levelStats(level).drillBestDist = drillBestDist;
    levelStats(level).lowerBound = lowerBound;
    
    
    %% BnB 7 - Visualizations (3 figures)
    
    if doPlots
        %% second figure: visualize important quadrilaterals on grid:
        if debugVis
            importantQuadsFig = figure(11); clf
            [gridXs,gridYs] = meshgrid(ImGridAxisPoints.Xpositions,ImGridAxisPoints.Ypositions);
            plot(gridXs,gridYs,'.k'); hold on;
            
            axis([bx,Bx,by,By]);
            
            hold on;
            ps = [];
            ps(end+1) = plot(bestCornersX([1:4 1]),bestCornersY([1:4 1]),'g','linewidth',2);
            ps(end+1) = plot(GT_targetCornerXs([1:4 1]),GT_targetCornerYs([1:4 1]),'r','linewidth',2);
            ps(end+1) = plot(nnCornersX([1:4 1]),nnCornersY([1:4 1]),'c','linewidth',2);
            ps(end+1) = plot(drillBestXs([1:4 1]),drillBestYs([1:4 1]),'m:','linewidth',2);
            ps(end+1) = plot(targetCorners.x([1:4 1]),targetCorners.y([1:4 1]),'k:','linewidth',1);
            %     ps(end+1) = plot(snappedGTxs([1:4 1]),snappedGTys([1:4 1]),'y','linewidth',2);
            legend(ps,{'best','GT','Nearest','Drill','image 2'})
            title(['Important quadrilaterals and grid, at level ' num2str(level)])
            axis equal; axis tight
        end
        %% second figure: improtant quadrilaterals on images and histogram of scores
        
        GPM_result_fig = figure(22); clf
        
        % visualize on images
        subplot 221; imshow(I1); title('I1');
        subplot 222; imshow(I2); title('I2');
        
        hold on;
        ps = [];
        ps(end+1) = plot(bestCornersX([1:4 1]),bestCornersY([1:4 1]),'g','linewidth',2);
        ps(end+1) = plot(GT_targetCornerXs([1:4 1]),GT_targetCornerYs([1:4 1]),'r','linewidth',2);
        ps(end+1) = plot(nnCornersX([1:4 1]),nnCornersY([1:4 1]),'c','linewidth',2);
        ps(end+1) = plot(drillBestXs([1:4 1]),drillBestYs([1:4 1]),'m:','linewidth',2);
        legend(ps,{'best','GT','Nearest','Drill'})
        
        % scores histogram
        h_Hist = subplot(2,2,3);
        height = 100*mean(percentileDists);
        hist(min(percentileDists,prctile(percentileDists,90)),30); hold on;
        title(sprintf('histogram of %d scores. #toExpand: %d. Delta: %.1f. Buffer: %.1f. '...
            ,numel(percentileDists),nnz(percentileDists<thresh),GridStepSize,guarantyFactor*GridStepSize));
        ps = [];
        ps(end+1) = plot([bestDist,bestDist],[0,height],'-g','linewidth',3);
        ps(end+1) = plot([drillBestDist,drillBestDist],[0,height],'-m','linewidth',3);
        ps(end+1) = plot([GT_bestDist,GT_bestDist],[0,height],':r','linewidth',3);
        ps(end+1) = plot([NN_bestDist,NN_bestDist],[0,height],'-c','linewidth',3);
        ps(end+1) = plot([thresh,thresh],[0,height],'-b','linewidth',3);
        ps(end+1) = plot([lowerBound,lowerBound],[0,height],'-k','linewidth',3);
        legend(ps,{['best ' num2str(bestDist,'%.1f')],['drill ' num2str(drillBestDist,'%.1f')],['GT ' num2str(GT_bestDist,'%.1f')],...
            ['Nearest ' num2str(NN_bestDist,'%.1f')],['thresh ' num2str(thresh,'%.1f')],['lower bound ' num2str(lowerBound,'%.1f')]})
        
        subplot 224; hold on;
        title('convergence')
        LS = levelStats(1:level);
        hs = plot(1:level,[[LS.thresh]',[LS.bestDist]',[LS.drillBestDist]',[LS.lowerBound]'],'-*');
        plotRng = [0.5 level+0.5];
        hs(end+1) = plot(plotRng,[GT_bestDist,GT_bestDist],':k');
        xlabel('B&B iteration');
        ylabel('error');
        legend(hs,{'upperBound','bestDist','drillBestDist','lowerBound','GT'})
        xlim(plotRng);
        set(gca,'XTick',1:level)
        
        
        %% third figure: IRE - Showing also number of expansions (for inlier rate determining)
        
        h_IRE = figure(33); clf
        subplot 131; plot(linspace(0,100,length(orig_GT_distances)),sort(orig_GT_distances)); hold on
        title(sprintf('"GT" transformation: \nsorted errors'),'fontweight','bold');
        plot(100*[chosenPercentile,chosenPercentile],[0,100],'r:','linewidth',3)
        xlabel('prcentile of scores')
        ylabel('error')
        ylim([0,300]);
        xlim([0,100]);
        
        subplot 132; hold on
        title(sprintf('bounds on: \n"best error per percentile"'),'fontweight','bold');
        plot((1:numMatches)/numMatches,bestVals ,'k','linewidth',2)
        plot((1:numMatches)/numMatches,bestVals+guarantyFactor*GridStepSize ,'r')
        plot((1:numMatches)/numMatches,max(0,bestVals-guarantyFactor*GridStepSize) ,'g')
        plot([chosenPercentile,chosenPercentile],[0,100],'r:','linewidth',3)
        xlabel('prcentile of scores')
        ylabel('error')
        %     ylim([0,max(bestVals)+guarantyFactor*GridStepSize])
        ylim([0,300]);
        xlim([0,1]);
        
        subplot 133; hold on
        title(sprintf('# of transformations at \na bounded distance from "best"'),'fontweight','bold');
        plot((1:numMatches)/numMatches,nExpansions)
        xlabel(['prcentile of scores. Chosen: ' num2str(chosenPercentile,'%.2f')])
        ylabel('# to be expanded')
        ylim([0,prctile(nExpansions,97)])
        xlim([0,1])
        plot([chosenPercentile,chosenPercentile],[0,prctile(nExpansions,90)/3],'r:','linewidth',3)
        
        saveas(h_IRE,[OutputDir 'IRE_round_' num2str(level) '_need_' num2str(needToEvaluateIRE) '.fig']);
        saveas(h_IRE,[OutputDir 'IRE_round_' num2str(level) '_need_' num2str(needToEvaluateIRE) '.png']);
    else
        GPM_result_fig = [];
    end
    
    %% BnB 8 - EXPANSION (incl. 1 figure)
    
    goodInds = find(percentileDists <= thresh);
    
    if isempty(goodInds)
        break
    end
    
    FilteredNetOfQuadrilaterals_indeces = NetOfQuadrilaterals_indices(:,goodInds);
    FilteredNetOfQuadrilaterals = NetOfQuadrilaterals(:,goodInds);
    
    % expand independently for factor 2 using unique
    
    factor = 2;
    neigh = 4; %8
    expansionRad = 3;
    if neigh == 8
        expansionRate = expansionRad^8;
    else % neigh == 4
        expansionRate = 5^4;
    end
    reppedIndicesMat = kron(factor*FilteredNetOfQuadrilaterals_indeces,ones(1,expansionRate,'int16'));
    perturbMat = GetPerturbMat('int16',expansionRad,neigh); % 8 x 9^4 or 8 x 5^4
    reppedPerturbMat = repmat(perturbMat,[1,length(goodInds)]);
    
    expandedIndices = reppedIndicesMat + reppedPerturbMat;
    NetOfQuadrilaterals_indices = unique(expandedIndices','rows')';
    
    GridStepSize = 0.5*GridStepSize;
    
    vecXs = [(ImGridAxisPoints.Xpositions - GridStepSize) ImGridAxisPoints.Xpositions(end) + GridStepSize];
    xpos = zeros(1,1+2*length(ImGridAxisPoints.Xpositions));
    xpos(1:2:end) = vecXs;
    xpos(2:2:end) = ImGridAxisPoints.Xpositions;
    ImGridAxisPoints.Xpositions = xpos;
    
    vecYs = [(ImGridAxisPoints.Ypositions - GridStepSize) ImGridAxisPoints.Ypositions(end) + GridStepSize];
    ypos = zeros(1,1+2*length(ImGridAxisPoints.Ypositions));
    ypos(1:2:end) = vecYs;
    ypos(2:2:end) = ImGridAxisPoints.Ypositions;
    ImGridAxisPoints.Ypositions = ypos;
    
    NetOfQuadrilaterals = zeros(size(NetOfQuadrilaterals_indices),'int16');
    NetOfQuadrilaterals(1:2:end,:) = ImGridAxisPoints.Xpositions(NetOfQuadrilaterals_indices(1:2:end,:));
    NetOfQuadrilaterals(2:2:end,:) = ImGridAxisPoints.Ypositions(NetOfQuadrilaterals_indices(2:2:end,:));
    
    
    % Expansion Figure: visualize the filtered grid-map and the new grid-map over the new grid:
    if debugVis
        newAndOldGridsFig = figure(44); clf; hold on;
        [gridXs_new,gridYs_new] = meshgrid(ImGridAxisPoints.Xpositions,ImGridAxisPoints.Ypositions);
        % plotting the new grid
        plot(gridXs_new,gridYs_new,'.k');
        % plotting the old grid
        plot(gridXs_new(2:2:end,2:2:end),gridYs_new(2:2:end,2:2:end),'.k','markersize',16);
        % setting axis according to dest. image size
        axis([bx,Bx,by,By]);
        for i = 1 : 4
            hhh = [];
            prev_unique_corner_coords = unique(FilteredNetOfQuadrilaterals(2*i-1:2*i,:)','rows')';
            hhh(end+1) = plot(prev_unique_corner_coords(1,:),prev_unique_corner_coords(2,:),'or','markersize',10);
            unique_corner_coords = unique(NetOfQuadrilaterals(2*i-1:2*i,:)','rows')';
            hhh(end+1) = plot(unique_corner_coords(1,:),unique_corner_coords(2,:),'Xb','markersize',10);
            hhh(end+1) = plot(GT_targetCornerXs(i),GT_targetCornerYs(i),'xr','markersize',10,'linewidth',3);
            hhh(end+1) = plot(nnCornersX(i),nnCornersY(i),'xc','markersize',10,'linewidth',3);
            hhh(end+1) = plot(snappedGTxs(i),snappedGTys(i),'xy','markersize',10,'linewidth',3);
            hhh(end+1) = plot(bestCornersX(i),bestCornersY(i),'xg','markersize',10,'linewidth',3);
            hhh(end+1) = plot(drillBestXs(i),drillBestYs(i),'xm','markersize',10,'linewidth',3);
            legend(hhh,{'good','expanded','GT','NN','snapped','best','drill'});
        end
        title(['Expansion process at level ' num2str(level)])
        axis equal; axis tight;
    end
    
    %% BnB 9 - Round summary and Stopping criteria
    
    fprintf('>>> Ended round %d: goodInds: %d, expanded: %d, unique: %d, GridStepSize: %.1f, innerChange: %.1f\n',...
        level,length(goodInds),size(expandedIndices,2),size(NetOfQuadrilaterals_indices,2),2*GridStepSize,2*GridStepSize*guarantyFactor)
    fprintf('>   drillBestDist: %.1f, bestDist: %.1f, GT_mean: %.1f, NN_mean: %.1f, thresh: %.1f, orig-thresh: %.1f, lowerBound: %.1f, \n',...
        drillBestDist,bestDist,GT_bestDist,NN_bestDist,thresh,origThresh,lowerBound)
    
    % criterion 1
    if drillBestDist < lowerBound + 0.1
        warning('drill error (%.2f) reached 0.1 from the lowerbound (%.2f)',drillBestDist,lowerBound);
        if ~goodIREevaluation
            fprintf('BREAKING BnB, by criterion 1 (drill score close to LB)\n');
            break
        end
    end
    
    % criterion 2
    if GridStepSize < 0.3
        fprintf('BREAKING BnB, by criterion 2 (small step size)\n');
        break
    end
    
    % criterion 3
    drillBestDists(end+1) = drillBestDist;
    if (level >= 4) && (drillBestDist > 0.98*drillBestDists(end-3)) && ~(drillBestDist > 1.02*drillBestDists(end-3))
        fprintf('BREAKING BnB, by criterion 3 (score hasn''t improved in 4 rounds)\n');
        break
    end
    
    % criterion 4
    if (level == maxBnBlevels)
        fprintf('BREAKING BnB, by criterion 4 (level == maxBnBlevels)\n');
        break
    end

    
end % levels of BnB

% evaluating H_drill on all points - check
% drillBestHomog = SquaresToHomography_mex(templateCorners.x',templateCorners.y',drillBestXs,drillBestYs);
% assert(isequal(drillBestHomog,H_drill));

distances_drill = ApplyHomographiesOnLocations_mex(H_drill,x1(1,:),x1(2,:),x2(1,:),x2(2,:)); % mex
if byPercentileAverage
    used_distances_drill = sort(distances_drill);
    used_distances_drill = bsxfun(@rdivide,cumsum(used_distances_drill,1),(1:length(x1))');
end

pStar = percentileToUse/robustPercentile;
nInliers = round(pStar*size(distances_drill,1)/100);
fprintf('---------------------------------------------------------------------------------------------------\n')
fprintf('GPM ran %d bnb levels, at error %.2f and LB %.2f, at robust percentile %.1f\n',...
    level,drillBestDist,lowerBound,percentileToUse);
fprintf('      Final percentile (pStar) is %.1f with %d inliers out of %d matches\n',...
    pStar,nInliers,numMatches);
fprintf('---------------------------------------------------------------------------------------------------\n')

status = 'success';

% wrap the result
GPMresult.status = status;
GPMresult.drillBestXs = drillBestXs;
GPMresult.drillBestYs = drillBestYs;
GPMresult.H_drill = H_drill;
GPMresult.distances_drill = distances_drill;
GPMresult.used_distances_drill = used_distances_drill;
GPMresult.used_GT_distances = GT_distances;
if exist('USAC_distances','var');
    GPMresult.used_USAC_distances = used_USAC_distances;
end
GPMresult.drillBestDist = drillBestDist;
GPMresult.pStarRobust = percentileToUse;
GPMresult.pStar = pStar; 
GPMresult.nInliers = nInliers;
GPMresult.lowerBound = lowerBound;
GPMresult.BnBlevels = level;
GPMresult.GPM_result_fig = GPM_result_fig;




