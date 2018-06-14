function [img_template, start_frame, last_frame, options, datapath] = bsmc_loadDataset(dataset)

% global mTracksAll;
% global mLabelAll;
basedir = strrep(fileparts(mfilename('fullpath')),'\','/');

options.lookahead = 5;
options.smoothSize = 5;
res = 3;


datapath = [ basedir, '/../../Data/moseg_dataset/', dataset];
img_template = [datapath, '/', dataset, '_%02d.ppm'];
%[start_frame, last_frame] = getFrameNums(datapath, 'jpg', [dataset, '_%02d.jpg']);

[start_frame, last_frame] = getFrameNums(datapath, 'ppm', [dataset,'_%02d.ppm']);

options.addNew = 1;
options.sigma_traj = 10;
options.description = 'Bg+Fg modelling, Window size = 7, window gaussian propability based on motion, post GC';
options.method = 'kde';
options.submethod = 'bg+fg';
options.kde_n = 10;
options.kde_thresh = 5e-7; % For people 1 sequence it seems we need to increase this threshold from 1e-9 to say 2e-9 or 1e-8 otherwise we loose most of the legs
options.kde_start_eval = 5;
options.win_size = 7;
options.colorspace = 'rgs';
options.postGC = 1;
options.borderInitNeighbor = 0;
options.label_prior = 0;
options.motion_window = 1;

% Example lookahed = 3 maximum smooth size 4
% or minimum lookahead is smooth size - 1
assert(options.smoothSize <= (options.lookahead+1));
end


function [start_frame, last_frame] = getFrameNums(datapath, extension, fileTemplate)
% Find the number of frame automaticaly
files = dir([datapath '/*.' extension]);
last_frame = 0;
start_frame = 1000000;
for i=1:length(files)
    frameNo = sscanf(files(i).name, fileTemplate);
    last_frame = max(last_frame, frameNo);
    start_frame = min(start_frame, frameNo);
end
end
