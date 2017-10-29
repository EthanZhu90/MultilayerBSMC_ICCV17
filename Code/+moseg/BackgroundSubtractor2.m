classdef BackgroundSubtractor2 < handle
    properties
        startframe;
        endframe;
        curframe;
        mosegframe;
      
        mosegdriver;
        Initialized = false;
        densesegmenter = [];
        Options;
    end
    
    methods
        function obj = BackgroundSubtractor2(varargin)

            if (mod(length(varargin),2) ~=0)
                error('Invalid number of arguments');
            end
            argused = true(1,length(varargin) ./ 2);            
            obj.startframe = 1;
            
            for i=1:2:length(varargin)
                switch(varargin{i})
                    case 'StartFrame'
                        obj.startframe = varargin{i+1};
                        argused((i+1)/2) = false;
                    case 'MosegAlgorithm'
                        varargin{i} = 'Algorithm';
                        argused((i+1)/2) = false;
                    otherwise
                        argused((i+1)/2) = false;
                end
            end
            obj.curframe = obj.startframe-1;
            
            mask = reshape([~argused; ~argused], 1, []);
            % Send remaining argumnets to MosegDriver2
            obj.mosegdriver = moseg.MosegDriver2(varargin{mask});
            obj.Initialized = false;
            
            %option for backgroundSubtractor
            options.smooth_lk = 0.01;
            options.prior_weight = 1;
            options.panelty_weight = 0.1;
            options.clustersize_thresh = 0.01; % rate of the cluster size to the whole image.
            obj.Options = options;
        end
        
        function initialize(obj, seq)
            % Perform itialization
            obj.Initialized = false;
            FramestoProcess = min(length(seq.Frames), 150); %% max is 150 Frames
            obj.endframe = seq.Frames(1) - 1 + FramestoProcess; %length(seq.Frames);
            obj.mosegframe = seq.Frames(1) - 1;
            obj.curframe = seq.Frames(1) - 1;  %% startFrame index
            
            outpath = sprintf('../Results/DataFiles/%s', seq.SeqName);
            mosegOutpath = sprintf('../Results/DataFiles/%s', seq.SeqName);
            % Create output directories if necessary
            if (~exist(outpath, 'dir')); mkdir(outpath); end
            if (~exist(mosegOutpath, 'dir')); mkdir(mosegOutpath); end
            for i = 1:4
                obj.mosegframe = obj.mosegframe + 1;
                mosegFname = [mosegOutpath '/mosegState-' num2str(obj.mosegframe) '.mat'];
                im = seq.readImage(obj.mosegframe);
                obj.mosegdriver.step(im);
                mosegState =  obj.mosegdriver.mosegmenter.CurrState;
                save(mosegFname, 'mosegState');
            end
        end
        
        
        function step(obj, seq)
            hm = vision.MarkerInserter;
            % Perform sparse segmentation
            obj.curframe = obj.curframe+1;
            
            % Save currState
            outpath = sprintf('../Results/DataFiles/%s', seq.SeqName);
            respath = sprintf('../Results/Results/%s', seq.SeqName );
            mosegOutpath = sprintf('../Results/DataFiles/%s', seq.SeqName);
            
            % Create output directories if necessary
            if (~exist(outpath, 'dir')); mkdir(outpath); end
            if (~exist(mosegOutpath, 'dir')); mkdir(mosegOutpath); end
            if (~exist(respath, 'dir')); mkdir(respath); end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Perform sparse segmentation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.mosegframe < obj.endframe
                obj.mosegframe = obj.mosegframe + 1;
                mosegFname = [mosegOutpath '/mosegState-' num2str(obj.mosegframe) '.mat'];
                im = seq.readImage(obj.mosegframe);
                obj.mosegdriver.step(im);
                mosegState =  obj.mosegdriver.mosegmenter.CurrState;
                save(mosegFname, 'mosegState');
                futureMosegState = mosegState;
            else
                mosegFname = [mosegOutpath '/mosegState-' num2str(obj.endframe) '.mat'];
                load(mosegFname, 'mosegState'); %% load mosegState
                futureMosegState = mosegState;
            end
            
            future3 = obj.curframe+2; 
            if future3 <= obj.endframe
                mosegFname = [mosegOutpath '/mosegState-' num2str(future3) '.mat'];
                load(mosegFname, 'mosegState'); %% load mosegState
                futureMosegState3 = mosegState; %% future 3 frame; 
            else
                mosegFname = [mosegOutpath '/mosegState-' num2str(obj.endframe) '.mat'];
                load(mosegFname, 'mosegState'); %% load mosegState
                futureMosegState3 = mosegState; 
            end
            mosegFname = [mosegOutpath '/mosegState-' num2str(obj.curframe) '.mat'];
            load(mosegFname, 'mosegState'); %% load mosegState
            %%% futureMosegState: 5; futureMosegState3: 3; mosegFname: curr

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Perform dense segmentation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            im = seq.readImage(obj.curframe);
            if ~isempty(mosegState.bgclust)
                bgclust = mosegState.bgclust;
            else
                bgclust = futureMosegState3.bgclust;
            end
            fprintf('======== DenseSeg Frame %03d\n', obj.curframe);
            if(obj.Initialized==false)
                obj.densesegmenter = moseg.LayeredModel(obj.curframe);
                obj.densesegmenter.step(im, mosegState, futureMosegState, futureMosegState3, bgclust);
                obj.Initialized = true;
            else
                obj.densesegmenter.step(im, mosegState, futureMosegState, futureMosegState3, bgclust);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% the segmentation is done so far
            %%% the rest code is to save result of current frame: 
            %%% 1. dense segmentation 2. trajectory  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            state =  obj.densesegmenter.CurrState;
            outim = state.frame;
            outim = double(outim);
            
            %%% figure 1: dense segmentation %%%
            %label_color = {'blue', 'red', 'green', 'light blue', 'yellow', 'pink', 'light blue', 'yellow', 'pink'};
            cmap = [1,0,0; 0,0,1; 0,1,0; 0,1,1; 1,1,0; 1,0,1; 0,0.5,0.5; 0.5,0.5,0; 0.5,0,0.5] *255;
            idx = find([state.layers.label] ~= state.bgClust);
            if(~isempty(idx))
                cnt = 1;
                for fglayer = idx
                    layer_mask = (state.lmask == state.layers(fglayer).label); 
                    if cnt > size(cmap, 1)
                        cnt = cnt - 1;
                    end
                    chnls = cmap(cnt,:);
                    for chnl_idx = 1:3 % only three channels.
                        if (chnls(chnl_idx) ~=0)
                            outim(:,:,chnl_idx) = double(outim(:,:,chnl_idx)) .* ~layer_mask  + layer_mask * chnls(chnl_idx);% *0.4
                        end
                    end
                    cnt = cnt + 1;
                end
            end
            layer_mask = (state.lmask == state.bgClust);
            if(~isempty(layer_mask))
                for i = 1:3
                    outim_Ch = outim(:,:,i); 
                    outim_Ch(layer_mask) =  outim_Ch(layer_mask) * 0.5; 
                    outim(:,:,i) = outim_Ch;
                end
            end
            outim = uint8(outim);
            outim = insertText(outim, [15 20],  num2str(obj.curframe), 'TextColor', 'white', 'BoxOpacity',0.0, 'FontSize',30);
            imwrite(outim, [respath '/MultiSeg-' num2str(obj.curframe) '.png'],'png');
            obj.densesegmenter.CurrState.MultiSegIm = outim; 
            
            %%% figure 2: Trajectory with diff label %%%
            outim = seq.readImage(obj.curframe);
            blankim = uint8(255 * ones(size(outim))); 
            
            cmap = [0,0,0; 1,0,0; 0,0,1; 0,1,0; 0,1,1; 1,1,0; 1,0,1; 0,0.5,0.5; 0.5,0.5,0; 0.5,0,0.5;] * 255;
            color_occupy = zeros(size(cmap,1),1);
            for i = 1:length(state.layers)
                points = round(state.layers(i).traj_store.sparse_points); 
                hm.Shape = 'Circle';
                hm.Fill = true; 
                hm.Size = 2;
                %hm.BorderColor = 'Custom';
                hm.FillColor = 'Custom';
                %hm.CustomBorderColor = cmap(state.layers(i).label,:);
                if(state.layers(i).label > size(cmap,1))
                    free_idx = find(color_occupy == 0);
                    if ~isempty(free_idx)
                        color_idx = free_idx(1);
                        color_occupy(color_idx) = 1;
                    else
                        color_idx = 10; % too many cluster, extra ones set to 10
                    end
                else
                    color_idx = state.layers(i).label;
                    color_occupy(color_idx) = 1;
                end
                hm.CustomFillColor = cmap(color_idx,:); 
                outim = step(hm, outim, int16(points)');
                blankim  = step(hm, blankim, int16(points)');
                release(hm);
            end
            imwrite(outim, [respath '/BPTracking_C&F-' num2str(obj.curframe) '.png'],'png');
            obj.densesegmenter.CurrState.trajmap = blankim; 
            obj.densesegmenter.CurrState.AllTrajIm = outim; 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% save the dense segmenter states
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            state =  obj.densesegmenter.CurrState;
            fname = [outpath '/state-' num2str(obj.curframe) '.mat'];
            save(fname, 'state');
           
        end
        
        function plot(obj, seq, iframe)
            
            statefile = sprintf('../Results/DataFiles/%s/state-%d.mat',  seq.SeqName,iframe);
            respath   = sprintf('../Results/Results/%s', seq.SeqName );
            if (~exist(respath, 'dir')); mkdir(respath); end
            load(statefile);   
           
            clf;
            [ncols, nrows, ~] = size(state.frame);
            nr = 2;
            nc = 3;
            resize_figure( nrows * nc, ncols * nr);
            
            %%%divides the current figure into an nr-by-nc grid and creates an axes for a subplot in the position specified by p
            axes_subplot(nr, nc, 1, 1);
            imshow(state.AllTrajIm);
            axis image off;
            text(5, 10, 'AllTraj', 'color','white','FontWeight','bold');
            
            axes_subplot(nr, nc, 1, 2);
            imshow(state.MultiSegIm);
            axis image off;
            text(5, 10, 'MultiSeg', 'color','white','FontWeight','bold');
            
            axes_subplot(nr, nc, 1, 3);
            imagesc(state.lmask);
            colormap default; 
            axis image off;
            text(5, 10, 'AfterMorph', 'color','white','FontWeight','bold');
            
            axes_subplot(nr, nc, 2, 1);
            imagesc(state.debug_info.InitSegIm);
            axis image off;
            text(5, 10, 'InitalGC', 'color','white','FontWeight','bold');
            
            axes_subplot(nr, nc, 2, 2);
            imagesc(~state.debug_info.app_invalid_LsPixel_Com);
            axis image off;
            text(5, 10, 'Valid_Red', 'color','white','FontWeight','bold');
            
            
            axes_subplot(nr, nc, 2, 3);
            imagesc(state.lmaskBfMorph);
            axis image off;
            text(5, 10, 'lmaskBfMorph', 'color','white','FontWeight','bold');

            drawnow;
            savefile = sprintf([respath '/Combo-%d.png'], iframe);
            saveas(gcf,savefile,'png') ; 
          
        end
    end

end