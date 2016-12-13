function [fileName1,fileName2,GT_homog,useprctile,scaleFact] = LoadImages(pthData,name)

switch name


    case 'graffiti_viewpoint_frames_1_4';
        fileName1 = fullfile(pthData,'miko\graffiti_viewpoint\img1.bmp');
        fileName2 = fullfile(pthData,'miko\graffiti_viewpoint\img4.bmp');
        GT_homog = dlmread(fullfile(pthData,'miko\graffiti_viewpoint\H1to4p'));
%         useprctile = 27;
        scaleFact = 2; % 1.4;

    case 'graffiti_viewpoint_frames_1_5';
        fileName1 = fullfile(pthData,'miko\graffiti_viewpoint\img1.bmp');
        fileName2 = fullfile(pthData,'miko\graffiti_viewpoint\img5.bmp');
        GT_homog = dlmread(fullfile(pthData,'miko\graffiti_viewpoint\H1to5p'));
%         useprctile = 5;
        scaleFact = 3; % 1.4;


    otherwise
        error('LoadImages: no such case!');

end