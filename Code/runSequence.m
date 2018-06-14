function [] = runSequence(segname)
%%% Path and Setup %%%%%%%%%%%%%%%%%%%
clc;
%clear; 

if(~exist('segname', 'var'))
    segname = 'cars1'; 
end

restoredefaultpath;
fprintf('Processing Sequence %s..........\n\n', segname); 
addpath(genpath('./BSMC/'));

ds = moseg.BroxDataset;
seq = ds.getSequence(segname);


detector = moseg.GridDetector(8);  %%stepsize is 8, 4 is too dense, 
flowpath = [seq.SeqRoot '/OpticalFlow/'];
tracker = moseg.LDOFTracker(flowpath);

algopt = moseg.MosegAlgLblPropOptions; %% store the option parameter
algopt.sigma_noise = 0.5;
alg = moseg.MosegAlgLblProp(1, 'Options', algopt);

bs = moseg.BackgroundSubtractor2('MosegAlgorithm', alg, 'Tracker', tracker, 'Detector', detector);

% % % process frames in online way
% % bs.initialize(seq);  
% % for i= seq.Frames
% %     bs.step(seq);
% %     % uncomment if you want check the intermidiate result
% %     %bs.plot(seq, i);
% %     %pause();
% % end

%%% generate visualized result, and score. 
gtpath  = sprintf('../Data/groundtruth/%s_gt.mat', segname); 
outpath = sprintf('../Results/DataFiles/%s', segname); 
respath = sprintf('../Results/Results/%s', segname );
mosegOutpath = sprintf('../Results/DataFiles/%s', segname);

[img_template, start_frame, last_frame, options] = bsmc_loadDataset(segname);
bsmc_exportResults(img_template, start_frame, last_frame, mosegOutpath, outpath, respath, gtpath, options, 'stats2'); 
bsmc_exportResults(img_template, start_frame, last_frame, mosegOutpath, outpath, respath, gtpath, options, 'evalMultiLabel');

