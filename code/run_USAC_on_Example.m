% USAC script

function resultUSAC = run_USAC_on_Example(examplePth,img1Filename,img2Filename)

%% 1] run
codeRoot = fileparts(mfilename('fullpath'));
usac_exe = fullfile(codeRoot,'third_party\USAC\msvc\USAC\Release\RunSingleTest.exe');
homogFlag = '1';
% cfgPath = 'D:\Dropbox\EpipolarMatching\data\mikolycic\graffiti_viewpoint\frames_1_5\example.cfg';
cfgPath = [examplePth 'example.cfg'];

cmd = [usac_exe ' ' homogFlag ' ' cfgPath];

[status,cmdout] = system(cmd,'-echo');

%% show results

%% paths and constants
%     clear;
orig_pts_file = 'orig_pts.txt';
inliers_file = 'inliers.txt';
homog_file = 'H.txt';

skip = 5;   % show subset of inliers based on skip size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read in images
im1 = imread(img1Filename);
im2 = imread(img2Filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read in original data points
fid_o = fopen(fullfile(examplePth, orig_pts_file), 'r');
num_pts = str2num(fgetl(fid_o));
m1 = zeros(2, num_pts);
m2 = zeros(2, num_pts);
for i = 1:num_pts
    temp = textscan(fgetl(fid_o), '%s');
    m1(1, i) = str2num(temp{1,1}{1});
    m1(2, i) = str2num(temp{1,1}{2});
    m2(1, i) = str2num(temp{1,1}{3});
    m2(2, i) = str2num(temp{1,1}{4});
end
fclose(fid_o);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read in inlier data
if ~exist(fullfile(examplePth, inliers_file),'file')
    inliers = zeros(length(m1),1);
    inliers_ind = [];
else
    inliers = textread(fullfile(examplePth, inliers_file));
%     if ~any(inliers)
%         inliers(:) = 1;
%     end
    inliers_ind = find(inliers > 0);
    n = randperm(length(inliers_ind));
    inliers_ind = inliers_ind(n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read in result homography
if exist(fullfile(examplePth, homog_file),'file')
    H_usac_T = textread(fullfile(examplePth, homog_file),'%f');
else
    fprintf('run_USAC_on_Example:NoHomogFile - homography file not found on "%s"',examplePth)
    H_usac_T = zeros(1,9);
end
H_usac = reshape(H_usac_T,3,3)';
if all(H_usac==0)
    H_usac = 0.5*eye(3);
    H_usac(9) = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% display images side by side and plot inliers
I = [];
padding = 15;
[M1,N1,K1]=size(im1);
[M2,N2,K2]=size(im2);
N3=N1+N2+padding;
M3=max(M1,M2);
oj=N1+padding;
oi=0;
I(1:M3, 1:N3, 1) = 255; I(1:M3, 1:N3, 2) = 255; I(1:M3, 1:N3, 3) = 255;
I(1:M1,1:N1,:) = im1;
I(oi+(1:M2),oj+(1:N2),:) = im2;

nMatches = length(inliers);
nInliers = length(inliers_ind);
prctInliers = 100*nInliers/nMatches;

h = figure(88);
title_str = {strrep(examplePth,'_',' ') ; sprintf('found %d inliers out of %d matches (%.1f%%)',nInliers,nMatches,prctInliers)};
imshow(uint8(I)); title(title_str); hold on; axis image; axis off;
for m = 1:skip+1:length(inliers_ind)
    plot([m1(1,inliers_ind(m)) m2(1,inliers_ind(m))+oj], [m1(2,inliers_ind(m)) m2(2,inliers_ind(m))+oi],'go-','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','g');
end
% saving
saveas(h,sprintf('%sUSAC_result__in_%d__m_%d__p_%.1f.png',examplePth,nInliers,nMatches,prctInliers));
   
resultUSAC.nMatches = nMatches;
resultUSAC.nInliers = nInliers;
resultUSAC.inlierInds = inliers_ind;
resultUSAC.prctInliers = prctInliers;
resultUSAC.status = status;
resultUSAC.cmdout = cmdout;
resultUSAC.H_usac = H_usac;

