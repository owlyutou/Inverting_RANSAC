function PerspectivePointMatching(runName,imageNames,datasetName,minimalSiftThreshold,...
    Noise,noiseSigmaInPixels,GPMparams,USACparams)

% close all
% dbstop if error


[~,pthData] = startup;
pthInit = pwd;

InlierEstimationFiguresDir = '../InlierEstimationFigures/';
mkdir(InlierEstimationFiguresDir);


% w = warning('query','last')
warning('off','images:imshow:magnificationMustBeFitForDockedFigure')
warning('off','MATLAB:iofun:UnsupportedEncoding')
warning('off','MATLAB:MKDIR:DirectoryExists')
warning('off', 'Coder:MATLAB:singularMatrix')

rng(1234)

doPlots = false;

% GPM related
if ~exist('runName','var'), runName = 'local'; end
if ~exist('GPMparams','var'), GPMparams.byPercentileAverage = 1; GPMparams.robustPercentile = 0.9; end

% SIFT related
if ~exist('minimalSiftThreshold','var'), minimalSiftThreshold = 1.3; end
if ~exist('Noise','var'), Noise = 0; end
if ~exist('noiseSigmaInPixels','var'), noiseSigmaInPixels = 0; end
if ~exist('datasetName','var'), datasetName = 'MIKO'; end

% USAC related
if ~exist('USACparams','var'), USACparams.USAC_thresh = 2.0; USACparams.USAC_LO = 1; end


noiseName = sprintf('noise_%d',noiseSigmaInPixels);

siftMethod = ['vlfeat_sifts_' num2str(minimalSiftThreshold,'%.2f')]; % 'blogs_sifts'; %
siftMethod = strrep(siftMethod,'.','_');


%% Load images

if ~exist('imageNames','var')
    imageNames = {...
        'graffiti_viewpoint_frames_1_2',...
        'graffiti_viewpoint_frames_1_3',...
        'graffiti_viewpoint_frames_1_4',... % 27%
        'graffiti_viewpoint_frames_1_5',... % 5% :-|
        'graffiti_viewpoint_frames_1_6',... % ?
        'wall_viewpoint_frames_1_2',...
        'wall_viewpoint_frames_1_3',... % 83%
        'wall_viewpoint_frames_1_4',... % 63%    'wall_viewpoint_frames_1_6'% only 4% inliers,...
        'wall_viewpoint_frames_1_5',... % 37%
        'wall_viewpoint_frames_1_6',... % 4% :-|
        'trees_blur_frames_1_2',...
        'trees_blur_frames_1_3',...
        'trees_blur_frames_1_4',...
        'trees_blur_frames_1_5',...
        'trees_blur_frames_1_6',... % 17%
        'ubc_compression_frames_1_6',... % 38%
        'bark_Z_R_frames_1_2',...
        'bark_Z_R_frames_1_3',...
        'bark_Z_R_frames_1_4',...
        'bark_Z_R_frames_1_5',... % 30%
        'bark_Z_R_frames_1_6',... % ?
        'bikes_blur_frames_1_6',...% ?
        'light_frames_1_6'% ?
        };
end

% imageNames = {'graffiti_viewpoint_frames_1_2'};
% imageNames = {'graffiti_5_viewpoint_frames_1_2'};

for imgInd = 1 : length(imageNames)
    
    close all
    
    %% 1] load the images and matches
    name = imageNames{imgInd};
    if isequal(name(1:4),'USAC')
        pthData = [pwd '\third_party\USAC\data\homog'];
    end
    %
    [img1Filename,img2Filename,GT_homog] = LoadImages(pthData,name);
    if isempty(GT_homog) || all(GT_homog(:)==0) % this happens a lot with USAC
        GT_homog = eye(3);
    end
    
    
    % get directories and scaleFact (different whether using USAC precalculated mathces or MIKO new ones)
    [exampleDir,OutputDir,namePrefix,nameSuffix,GPMparams.scaleFact] = ...
        GetDirectoriesAndScaleFact(name,runName,siftMethod,noiseName,img1Filename);
    
    % compute/load matches and
    [I1,I2,x1,x2,ordering,numMatches] = ComputeOrLoadMatches(...
        siftMethod,minimalSiftThreshold,Noise,noiseSigmaInPixels,img1Filename,img2Filename,OutputDir);
    
    
    %% 2] Run USAC
    
    % run USAC
    resultUSAC = RunUSAC(datasetName,noiseName,namePrefix,nameSuffix,...
        siftMethod,exampleDir,x1,x2,ordering,img1Filename,img2Filename,USACparams);
    
    
    
    %% 3] Evaluate the USAC and GT results
    
    % image sizes
    [h1,w1,d1] = size(I1);
    [h2,w2,d2] = size(I2);
    
    % template corners
    templateCorners.x = [1,1,w1,w1];
    templateCorners.y = [1,h1,h1,1];
    targetCorners.x = [1,1,w2,w2];
    targetCorners.y = [1,h2,h2,1];
    
    % evaluate USAC
    USAC_homog = resultUSAC.H_usac;
    USAC_distances = ApplyHomographiesOnLocations_mex(reshape(USAC_homog',9,1),x1(1,:),x1(2,:),x2(1,:),x2(2,:)); % mex
    USACpercentile = resultUSAC.prctInliers;
    zVec = [-1,-1,-1,-1];
    [~,USAC_targetCornerXs,USAC_targetCornerYs] = ...
        ApplyHomographiesOnLocationsFull_mex(reshape(USAC_homog',9,1),templateCorners.x,templateCorners.y,zVec,zVec); % mex
    
    % evaluate GT
    orig_GT_distances = ApplyHomographiesOnLocations_mex(reshape(GT_homog',9,1),x1(1,:),x1(2,:),x2(1,:),x2(2,:)); % mex
    zVec = [-1,-1,-1,-1];
    [~,GT_targetCornerXs,GT_targetCornerYs] = ...
        ApplyHomographiesOnLocationsFull_mex(reshape(GT_homog',9,1),templateCorners.x,templateCorners.y,zVec,zVec); % mex
    
    
    if doPlots
        GTvsUSACpdfs_Fig = figure; hold on
        hUS = cdfplot(USAC_distances);
        hGT = cdfplot(orig_GT_distances);
        set(hGT,'color','r');
        legend([hUS,hGT],'USAC','GT');
    end
    
    %% 4] Run GPM
    
    %%%%%%
    if exist([OutputDir 'GPM_results.mat'],'file')
        load([OutputDir 'GPM_results.mat'],'GPMresult','GT_targetCornerXs','GT_targetCornerYs',...
            'drillBestXs','drillBestYs','drillBestDist','lowerBound',...
            'H_drill','GT_homog','pStar','BnBlevels',...
            'sDistances_drill','sDrillInds','nInliers','inlierIndsDrill');
    end
   
    
    if exist([OutputDir 'GPM_results.mat'],'file')
        continue % this is a sign that computation is done, so continue to next pair
    end
    
    GPMresult = RunGPM(GPMparams,x1,x2,OutputDir,...
        GT_targetCornerXs,GT_targetCornerYs,targetCorners,templateCorners,orig_GT_distances,I1,I2,USAC_distances);
    
    
    if ~strcmp(GPMresult.status,'success')
        disp(GPMresult.status);
        % save some of our results in a mat-file
        save([OutputDir 'GPM_results.mat'],'GPMresult','GT_targetCornerXs','GT_targetCornerYs','GT_homog');
        continue
    end
    
    % unwrap the result
    drillBestXs = GPMresult.drillBestXs;
    drillBestYs = GPMresult.drillBestYs;
    H_drill = GPMresult.H_drill; %#ok<NASGU>
    distances_drill = GPMresult.distances_drill;
    used_distances_drill = GPMresult.used_distances_drill;
    used_GT_distances = GPMresult.used_GT_distances;
    used_USAC_distances = GPMresult.used_USAC_distances;
    drillBestDist = GPMresult.drillBestDist;
    pStar = GPMresult.pStar;
    pStarRobust = GPMresult.pStarRobust;
    nInliers = GPMresult.nInliers;
    lowerBound = GPMresult.lowerBound; %#ok<NASGU>
    BnBlevels = GPMresult.BnBlevels; %#ok<NASGU>
    GPM_result_fig = GPMresult.GPM_result_fig;
    
    
    %% 5] Post processing - (mainly: 3 figures)
    
    % sorting all distances
    [sDistances_drill,sDrillInds] = sort(distances_drill);
    [sGT_distances,sGTInds] = sort(orig_GT_distances);
    [sUSAC_distances,sUSACInds] = sort(USAC_distances);
    
    inlierIndsDrill = sDrillInds(1:nInliers);
    
    %% Fig 0: matches on images (each of GPM and USAC, compared to GT)
    if doPlots
        N2show = 200;
        outlierIndsDrill = setdiff((1:length(x1))',inlierIndsDrill);
        fact = max(length(inlierIndsDrill),length(outlierIndsDrill))/N2show;
        if fact < 1, fact = 1; end
        inlierInds = inlierIndsDrill(randsample(length(inlierIndsDrill),ceil(length(inlierIndsDrill)/fact)));
        outlierInds = outlierIndsDrill(randsample(length(outlierIndsDrill),ceil(length(outlierIndsDrill)/fact)));
        
        % subplots 1-2: GPM
        allResultsOnImg_1 = figure; imshow(I1); hold on;
        plot(x1(1,inlierInds),x1(2,inlierInds),'.b','markersize',15);
        plot(x1(1,outlierInds),x1(2,outlierInds),'.r','markersize',15);
        % title(['I1: "GPM" inliers(' num2str(pStar,'%.1f') '% - BLUE) outliers (RED), with error: ' num2str(drillBestDist,'%.1f')]);
        legend({'Inliers','Outliers'});
        allResultsOnImg_2 = figure; imshow(I2); hold on;
        plot(x2(1,outlierInds),x2(2,outlierInds),'.r','markersize',15);
        plot(x2(1,inlierInds),x2(2,inlierInds),'.b','markersize',15);
        ps = [];
        ps(end+1) = plot(GT_targetCornerXs([1:4 1]),GT_targetCornerYs([1:4 1]),'g','linewidth',4);
        ps(end+1) = plot(drillBestXs([1:4 1]),drillBestYs([1:4 1]),'m','linewidth',4);
        % ps(end+1) = plot(USAC_targetCornerXs([1:4 1]),USAC_targetCornerYs([1:4 1]),'r','linewidth',2);
        legend(ps,{'GT','GMD'});
        % legend(ps,{'GT','GPM','USAC'});
        % title(['I2: "GPM" inliers(' num2str(pStar,'%.1f') '% - BLUE) outliers (RED), with error: ' num2str(drillBestDist,'%.1f')]);
        %
        
        saveas(allResultsOnImg_1,[OutputDir 'results_on_imgs_1.eps'],'epsc');
        saveas(allResultsOnImg_1,[OutputDir 'results_on_imgs_1.fig']);
        saveas(allResultsOnImg_1,[OutputDir 'results_on_imgs_1.png']);
        saveas(allResultsOnImg_2,[OutputDir 'results_on_imgs_2.eps'],'epsc');
        saveas(allResultsOnImg_2,[OutputDir 'results_on_imgs_2.fig']);
        saveas(allResultsOnImg_2,[OutputDir 'results_on_imgs_2.png']);
    end

    
    
    %% Fig 1: matches on images (each of GPM and USAC, compared to GT)
    
    allResultsOnImg = figure;
    % subplots 1-2: GPM
    subplot 221; imshow(I1); hold on;
    plot(x1(1,:),x1(2,:),'.r');
    plot(x1(1,inlierIndsDrill),x1(2,inlierIndsDrill),'.b');
    title(['I1: "GPM" inliers(' num2str(pStar,'%.1f') '% - BLUE) outliers (RED), with error: ' num2str(drillBestDist,'%.1f')]);
    subplot 222; imshow(I2); hold on;
    plot(x2(1,:),x2(2,:),'.r');
    plot(x2(1,inlierIndsDrill),x2(2,inlierIndsDrill),'.b');
    ps = [];
    ps(end+1) = plot(GT_targetCornerXs([1:4 1]),GT_targetCornerYs([1:4 1]),'g','linewidth',2);
    ps(end+1) = plot(drillBestXs([1:4 1]),drillBestYs([1:4 1]),'b','linewidth',2);
    ps(end+1) = plot(USAC_targetCornerXs([1:4 1]),USAC_targetCornerYs([1:4 1]),'r','linewidth',2);
    legend(ps,{'GT','GPM','USAC'});
    title(['I2: "GPM" inliers(' num2str(pStar,'%.1f') '% - BLUE) outliers (RED), with error: ' num2str(drillBestDist,'%.1f')]);
    %
    
    
    % subplots 1-2: USAC
    inlierIndsUSAC = resultUSAC.inlierInds;
    subplot 223; imshow(I1); hold on;
    plot(x1(1,:),x1(2,:),'.r');
    plot(x1(1,inlierIndsUSAC),x1(2,inlierIndsUSAC),'.b');
    title(['I1: "GT" inliers(' num2str(USACpercentile,'%.1f') '% - BLUE) outliers (RED), with error: 2']);
    subplot 224; imshow(I2); hold on;
    plot(x2(1,:),x2(2,:),'.r');
    plot(x2(1,inlierIndsUSAC),x2(2,inlierIndsUSAC),'.b');
    ps = [];
    ps(end+1) = plot(GT_targetCornerXs([1:4 1]),GT_targetCornerYs([1:4 1]),'g','linewidth',2);
    ps(end+1) = plot(drillBestXs([1:4 1]),drillBestYs([1:4 1]),'b','linewidth',2);
    ps(end+1) = plot(USAC_targetCornerXs([1:4 1]),USAC_targetCornerYs([1:4 1]),'r','linewidth',2);
    legend(ps,{'GT','GPM','USAC'});
    title(['I2: "USAC" inliers(' num2str(USACpercentile,'%.1f') '% - BLUE) outliers (RED), with error: 2']);
    
    % draw different inliers:
    diffInds = setdiff(inlierIndsDrill,inlierIndsUSAC);
    subplot 221; plot(x1(1,diffInds),x1(2,diffInds),'oc','markersize',12,'linewidth',2);
    subplot 222; plot(x2(1,diffInds),x2(2,diffInds),'oc','markersize',12,'linewidth',2);
    diffInds = setdiff(inlierIndsUSAC,inlierIndsDrill);
    subplot 223; plot(x1(1,diffInds),x1(2,diffInds),'oc','markersize',12,'linewidth',2);
    subplot 224; plot(x2(1,diffInds),x2(2,diffInds),'oc','markersize',12,'linewidth',2);
    
    %% Fig 2: histogram over our pecent of inliers
    if doPlots        
        figure; hist([sGT_distances(1:nInliers) sDistances_drill(1:nInliers) sUSAC_distances(1:nInliers)],100);
        legend('GT distances','GPM distances','USAC distances');
    end
    %% Fig 3: cdf of distances
    h_cdfs = figure; hold on
    hGT = cdfplot(sGT_distances); set(hGT,'color','g','linewidth',3,'displayname','GT');
    hGPM = cdfplot(sDistances_drill); set(hGPM,'color','b','linewidth',3,'displayname','GPM');
    hUSC = cdfplot(sUSAC_distances); set(hUSC,'color','r','linewidth',3,'displayname','USAC');
%     hUsedGpmDists = plot(used_distances_drill',(1/length(x1):1/length(x1):1),'b:','linewidth',2,'displayname','GPM used dists');
%     hUsedGtDists = plot(used_GT_distances',(1/length(x1):1/length(x1):1),'g:','linewidth',2,'displayname','GT used dists');
%     hUsedUsacDists = plot(used_USAC_distances',(1/length(x1):1/length(x1):1),'r:','linewidth',2,'displayname','USAC used dists');
    pGPM = plot([0,max(sGT_distances)],[0.01*pStar,0.01*pStar],'--k','linewidth',3,'displayname','GPM Inlier rate');
%     pGPM_rob = plot([0,max(sGT_distances)],[0.01*pStarRobust,0.01*pStarRobust],'--b','linewidth',2,'displayname','GPM robust inlier rate');
%     errorGPM = plot([drillBestDist,drillBestDist],[0,1],':b','linewidth',3,'displayname','GPM lowest error');
    legend('show','location','east');
    xlim([-0.5,min(40,drillBestDist*8)]);
    ylim([-0.05,1]);
    title(sprintf('GPM - robust error (%.2f) , robust pStar (%.1f) , pStar (%.1f)',drillBestDist,pStarRobust,pStar))
    
    
    %% Saving figures and workspace in an m-file
    
    if ishandle(GPM_result_fig), saveas(GPM_result_fig,[OutputDir 'GPM_result_fig.png']);end
    if ishandle(allResultsOnImg),saveas(allResultsOnImg,[OutputDir 'results_on_imgs.fig']); end
    if ishandle(h_cdfs),saveas(h_cdfs,[OutputDir 'cdf_comparisons.fig']);end
    % save some of our results in a mat-file
    save([OutputDir 'GPM_results.mat'],'GPMresult','GT_targetCornerXs','GT_targetCornerYs',...
        'drillBestXs','drillBestYs','drillBestDist','lowerBound',...
        'H_drill','GT_homog','pStar','BnBlevels',...
        'sDistances_drill','sDrillInds','nInliers','inlierIndsDrill');
    
    
end % loop over image pairs

%%
return


