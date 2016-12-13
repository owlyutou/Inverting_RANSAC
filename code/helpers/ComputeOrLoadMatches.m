% compute/load matches and scaleFact
function [I1,I2,x1,x2,ordering,numMatches] = ...
    ComputeOrLoadMatches(siftMethod,minimalSiftThreshold,Noise,noiseSigmaInPixels,img1Filename,img2Filename,OutputDir,pthBLOGS,pthInit)


%% compute or load sift matches by Blogs or VlFeat

siftFilename = [OutputDir siftMethod '.mat'];
if ~exist(siftFilename,'file')
    
    switch siftMethod(1:6)
        case 'blogs_' %%%%%%%%%%
            if ~strcmp(pwd,pthBLOGS),cd(pthBLOGS),end
            
            imwrite(rgb2gray(imread(img1Filename)),'BLOGS_onZuBuD_temp1.png')
            imwrite(rgb2gray(imread(img2Filename)),'BLOGS_onZuBuD_temp2.png')
            
            % sifts
            [I1, MD1, ML1, row1, col1] = sift('BLOGS_onZuBuD_temp1.png');  %Your first image in pgm format here
            [I2, MD2, ML2, row2, col2] = sift('BLOGS_onZuBuD_temp2.png');  %Your second image in pgm format here
            
            % matching
            [x1, x2, p_samp] = sift_match( MD1, ML1, MD2, ML2);
            numMatches = size(x1,2) ;
            ordering = zeros(numMatches,1);
            
            if ~strcmp(pthInit,pwd),cd(pthInit),end
            
        case 'vlfeat' %%%%%%%%%
            if exist('vl_sift','file')~=3
                error('please download and run vl_feat')
            end
            % read images
            I1 = imread(img1Filename); I2 = imread(img2Filename);
            
            % make single
            im1 = im2single(I1); im2 = im2single(I2);
            
            % make grayscale
            if size(im1,3) > 1, im1g = rgb2gray(im1) ; else im1g = im1 ; end
            if size(im2,3) > 1, im2g = rgb2gray(im2) ; else im2g = im2 ; end
            
            maxNumMatches = 500;
            
            % compute
            [f1,d1] = vl_sift(im1g);
            [f2,d2] = vl_sift(im2g);
            siftThreshold = 7;
            scores = [];
            while(numel(scores)<maxNumMatches && siftThreshold > minimalSiftThreshold)
                [matches, scores] = vl_ubcmatch(d1,d2,siftThreshold);
                siftThreshold = 0.95*siftThreshold;
                fprintf('siftThreshold: %.1f, numMatches: %d\n',siftThreshold,numel(scores));
            end
            
            
            % matches
            x1 = f1(1:2,matches(1,:)) ; x1(3,:) = 1 ;
            x2 = f2(1:2,matches(2,:)) ; x2(3,:) = 1 ;
            
            %% 3] take a random sample of matches and perhaps add some noise
            
            orig_x1 = x1;
            orig_x2 = x2;
            
            if length(x1) > maxNumMatches
                %                     [sortedScores,ordering] = sort(scores,'descend');
                %                     RS = ordering(1:maxNumMatches);
                RS = randsample(length(x1),maxNumMatches);
                x1 = x1(:,RS); x2 = x2(:,RS);
                scores = scores(RS);
            end
            numMatches = size(x1,2);
            
            if Noise
                x1(1:2,:) = x1(1:2,:) + noiseSigmaInPixels*randn(2,numMatches);
                x2(1:2,:) = x2(1:2,:) + noiseSigmaInPixels*randn(2,numMatches);
            end
            
            % ordering
            [sortedScores,ordering] = sort(scores,'ascend');
            ordering = ordering - 1; % zero-based
            
            
            
        otherwise
            error('no such sift method')
    end
    % save
    save(siftFilename,'I1','I2','x1','x2','ordering','numMatches','orig_x1','orig_x2','scores')
    
else
    % load
    load(siftFilename,'I1','I2','x1','x2','ordering','numMatches')
end
