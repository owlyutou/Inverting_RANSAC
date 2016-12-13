%%
clc
clear
close all
startup;

% 3 example image pairs
imageNames = {'graffiti_viewpoint_frames_1_4','graffiti_viewpoint_frames_1_5'};

Noise = 0;
noiseSigmaInPixels = 0;

GPMparams.byPercentileAverage = 1;
minimalSiftThreshold = 1.24;
USACparams.USAC_thresh = 1:31; % 2.0;  %
USACparams.USAC_LO = 1;

runName = 'MWE';
datasetName = 'Miko';

PerspectivePointMatching(runName,imageNames,datasetName,minimalSiftThreshold,Noise,noiseSigmaInPixels,GPMparams,USACparams);
