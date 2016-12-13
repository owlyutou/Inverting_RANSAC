% run USAC
function resultUSAC = RunUSAC(datasetName,noiseName,namePrefix,nameSuffix,...
    siftMethod,exampleDir,x1,x2,ordering,img1Filename,img2Filename,USACparams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unpacking params
USAC_thresh = USACparams.USAC_thresh;  %
USAC_LO = USACparams.USAC_LO;

%%
% USAC_thresh = 1:40;
for thresh = USAC_thresh
    fprintf('Running USAC on %s %s with threshold %d\n',namePrefix,nameSuffix,thresh);
    baseFileDir = [exampleDir 'USAC/'];
    threshSubdir = ['thresh_' num2str(thresh,'%d') '_LO_' num2str(USAC_LO,'%d')];
    matchfileDir = [baseFileDir threshSubdir '/'];
    resultUSACfile = [matchfileDir 'resultUSAC.mat'];
    if ~exist(resultUSACfile,'file')
        mkdir(matchfileDir);
        
        % 1] sorting-file
        sortingfileName = [matchfileDir 'sorting.txt'];
        
        fid = fopen(sortingfileName,'w');
        fprintf(fid,sprintf('%d\n',int16(ordering)));
        fclose(fid);
        
        % 2] match-file
        matchfileName = [matchfileDir 'orig_pts.txt'];
        
        fid = fopen(matchfileName,'w');
        nMatches = int16(size(x1,2));
        fprintf(fid,sprintf('%d\n',nMatches));
        allMatches = [x1(1:2,:)' x2(1:2,:)'];
        fprintf(fid,sprintf('%.2f %.2f %.2f %.2f\n',allMatches'));
        fclose(fid);
        
        % 3] config-file
        newCfgFile = [matchfileDir 'example.cfg'];
        codeRoot = strrep(fileparts(mfilename('fullpath')),'\','\\\\');
        codeRoot = strrep(codeRoot,'\\\\code','\\\\');
        
        copyfile('third_party\USAC\data\homog\template_example.cfg',newCfgFile);
        s1 = replaceinfile('D:\\\\Dropbox\\\\GPM\\\\code\\\\third_party\\\\USAC\\\\data\\\\homog\\\\test1\\\\',...
            [codeRoot 'results\\\\' datasetName '\\\\' siftMethod '\\\\' noiseName '\\\\' namePrefix '\\\\' nameSuffix '\\\\USAC\\\\' threshSubdir '\\\\'],newCfgFile );
        s2 = replaceinfile('inlier_threshold = 2.0;',['inlier_threshold = ' num2str(thresh,'%.1f') ';'],newCfgFile );
        
        if ~isempty(ordering) && ~all(ordering==0) % USAC can use PROSAC
            replaceinfile('\"UNIFORM\"','\"PROSAC\"',newCfgFile );
        end
        if ~USAC_LO
            replaceinfile('\"LOSAC\"','\"NO_LO\"',newCfgFile );
        end
        
        % run USAC
        resultUSAC = run_USAC_on_Example(matchfileDir,img1Filename,img2Filename);
        
        save(resultUSACfile,'resultUSAC');
    else
        load(resultUSACfile,'resultUSAC');
    end
end