classdef LayeredModel < handle
    properties
        Options;
        PrevState;
        CurrState;
        StartFrame;
    end
    
    methods
        function obj = LayeredModel(start_frame)
            obj.CurrState.frameNo = start_frame-1;
            obj.CurrState.layers = struct('label',{}, 'traj_store',{}, 'motion_model',{},'app_model', {}, 'prob_model',{}, 'mask',{}); 
            
            options.lrgNum = 100000000000000;
            options.init_app_gabp = 0;
            options.motion_subpixel = 0;
            
            % These affect the inialization of the appearance model.
            options.init_app_sigma_unary = 4;
            options.init_app_sigma_pairwise = 40;
            
            options.lookahead = 5;
            
            % Set defaults
            options.numObjects = 2;
            
            % Options for tracking.
            options.dense = 1;
            
            % Options for motion model inference.
            options.sig_edge = 0.5;
            options.sig_motion = 1;
            
            % Options for dense segmentation.
            options.sigma_app = 10;
            options.lambda = 66;
            options.boundary_pixels = 0;
            options.lambda_gc = 1;
            options.numMixtures = 3;
            
            % Options for appearance model prediction and update.
            % Amount of noise in color information between consecucative frames.
            % If low will depend more on the observation.
            options.sig_app2 = 10^2;
            
            % Options for clustering.
            options.dense_seg = 1;
            options.embed = 1;
            options.smooth_img = 0;
            options.nonrigid_thresh = 1;
            options.fg_selection_method = 'coordinates';
            options.pickFGOnce = 1;
            
            options.k_kmeans = 3;
            options.lookahead = 6;
            options.sigma_sp = 400;
            options.lambda_gc = 0.1;
            
            options.description = 'Bg+Fg modelling, Window size = 7, window gaussian propability based on motion, post GC=2';
            options.method = 'kde';
            options.submethod = 'bg+fg';
            options.kde_n = 10;
            options.kde_thresh = 5e-7;
            options.kde_start_eval = 5; %3; % determine the bginvalid and fginvalid 
            options.win_size = 5;  %7 win_size to compute post probability, replace the original 7 to 5 to speed up.
            options.colorspace = 'rgs';
            options.postGC = 2;    %1 is label as data cost, 2 is prob as data cost.  
            options.borderInitNeighbor = 0;
            options.label_prior = 0;
            options.motion_window = 1;
            
            options.minkdesig = 5;
            options.kde_shift = options.kde_thresh / 1.0001;
            options.trivial_cluster_thresh = 0.008; 
            %options.kde_n = 5;
            options.minTraj_rate = 0.005;
            options.SmallCompThresh = 50; 
            
           %%%%%%%%%% for graph_cut parameter%%%%%%%
            options.InGC_datacost_lambda=2;
            options.InGC_pairwise_lambda=8;
            options.InGC_labelcost_lambda=0.1; 
            
            options.datacost_lambda=2;
            options.pairwise_lambda=8;
            options.labelcost_lambda=0.1;
            
            obj.Options = options;
        end
        
        function step(obj, frame, mosegState, futureMosegState,futureMosegState3,  bgClust)
    
            label_color = {'blue', 'red', 'green', 'light blue', 'yellow', 'pink', 'light blue', 'yellow', 'pink'};
            obj.PrevState = obj.CurrState;
            
            obj.CurrState = struct;
            obj.CurrState.frameNo =  obj.PrevState.frameNo +1; 
            obj.CurrState.frame = frame;
            obj.CurrState.obs = bsmc_RGB2rgs(frame);
            obj.CurrState.bgClust = bgClust; 
            obj.CurrState.layers = struct('label',{},  'nframe', {}, 'traj_store',{}, 'motion_model',{},'app_model', {}, 'prob_model',{}, 'mask',{}); 
          
            if (isempty(mosegState.points))
                return;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% prepare the labeled trajectory 
            %%% combine the future 5 and 3 together. 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            future_memTrajIds5 = futureMosegState.trajIds(futureMosegState.lbls>0); %futureMosegState.memTrajIds; %
            future_lbls5 = futureMosegState.lbls(futureMosegState.lbls>0); %(futureMosegState.clust.lbls)'; %
           
            future_memTrajIds3 = futureMosegState3.trajIds(futureMosegState3.lbls>0); %futureMosegState.memTrajIds; %
            future_lbls3 = futureMosegState3.lbls(futureMosegState3.lbls>0); %(futureMosegState.clust.lbls)'; %
            
            [~, ia, ib] = intersect(future_memTrajIds5, future_memTrajIds3);
            future_memTrajIds = future_memTrajIds5(ia); 
            future_lbls = future_lbls5(ia); 
           
            [~, diff] = setdiff(future_memTrajIds5, future_memTrajIds3);
            future_memTrajIds = [future_memTrajIds, future_memTrajIds5(diff)]; 
            future_lbls = [future_lbls, future_lbls5(diff)]; 
            
            [~, diff] = setdiff(future_memTrajIds3, future_memTrajIds5);
            future_memTrajIds = [future_memTrajIds, future_memTrajIds3(diff)]; 
            future_lbls = [future_lbls, future_lbls3(diff)];
            
            % back propagation
            % future
            [~, ia, ib] = intersect(mosegState.trajIds, future_memTrajIds);
            lbls = future_lbls(ib);

            % current
            [~, diff] = setdiff(mosegState.trajIds, future_memTrajIds);
            albls = mosegState.lbls(diff); 
            c_lbls = albls(albls > 0); 
        
            c_points = mosegState.points(1:2, diff);
            obj.CurrState.Unlb_points =  c_points(1:2, albls == 0); 
          
            com_lbls = [lbls, c_lbls]; % labels storage 
            num_lbls = length(com_lbls); % total number of valid labels
            max_lbls = max(com_lbls); % the maximum index of labels
 
            valid_lbls = false(1,max_lbls);
            
            for i = 1: max_lbls
                if(sum(com_lbls == i)> obj.Options.minTraj_rate * num_lbls)
                    valid_lbls(i) = true; 
                end
            end
     
            if(valid_lbls(bgClust) == false)
                if sum(com_lbls == bgClust) > 0
                    fprintf('\nThe background does not contains enough trajectory! Force segmetation!...\n\n'); 
                    valid_lbls(bgClust) = true;
                else
                    fprintf('\nThere is no background cluster...\n\n'); 
                end
            end
            
            temp_lbls = 1: max_lbls;
            color_occupy = valid_lbls;
            if(length(color_occupy) > length(label_color))
                color_occupy = color_occupy(1:length(label_color));
            end
            Cnt = 0; 
            for i = temp_lbls(valid_lbls)
                
                [~, diff] = setdiff(mosegState.trajIds, future_memTrajIds);
                albls = mosegState.lbls(diff); 
                lbls = albls(albls > 0); 
                c_points = mosegState.points(1:2, diff);
                points = c_points(1:2, albls > 0); 
           
                current_sparse_points = round(points(1:2, (obj.matches_any(lbls,i))));
                % back propagation
                [~, ia, ib] = intersect(mosegState.trajIds, future_memTrajIds);
                lbls = future_lbls(ib);
                % Propagate labels back from future mosegState to current state.  
                sparse_points = round(mosegState.points(1:2, ia(obj.matches_any(lbls,i)))); 
                % combine the label from back propagation and the current label 
                sparse_points = [sparse_points, current_sparse_points]; 
               
                %%%%%%%%%%%%%%% prepare the m points
                [~, diff] = setdiff(mosegState.trajIds, future_memTrajIds);
                albls = mosegState.lbls(diff); 
                alen = mosegState.len(diff); 
                points = mosegState.points(1:4, diff);
                points = points(1:4, (albls > 0)&(alen >= 2)); 
                lbls = albls((albls > 0)&(alen >= 2)); 
              
                current_sparse_points_m = round(points(1:4, (obj.matches_any(lbls,i)))); % the point is same as sparse_points

                ib = ib(mosegState.len(ia) >=2);
                ia = ia(mosegState.len(ia) >=2);
                lbls = future_lbls(ib);             
                sparse_points_m = round(mosegState.points(1:4, ia(obj.matches_any(lbls,i))));    
                sparse_points_m = [sparse_points_m, current_sparse_points_m]; 
               
                Cnt = Cnt+1; 
                traj_store = struct('sparse_points', sparse_points, 'sparse_points_m', sparse_points_m); 
                
                %%% find the nframe
                idx = find([obj.PrevState.layers.label] == i); % find the corresponding layer
                if(isempty(idx))
                    obj.CurrState.layers(Cnt) =  struct('label', i, 'nframe', 1, 'traj_store', traj_store, 'motion_model', -1,'app_model', -1, 'prob_model',-1,  'mask', -1); 
                else
                    preNframe = obj.PrevState.layers(idx).nframe; 
                    obj.CurrState.layers(Cnt) =  struct('label', i, 'nframe', preNframe+1, 'traj_store', traj_store, 'motion_model', -1,'app_model', -1, 'prob_model',-1,  'mask', -1); 
                end
                
                %%% The IF statement below is Just for display, won't affect the segmentation result
                if(i > length(label_color))
                    free_idx = find(color_occupy==0);
                    if isempty(free_idx)
                        fprintf('Too many regions, no free color option!\n');
                        fprintf('The to-do-segmentation cluster with label %d and color %s   Traj #: %d \n', i, label_color{1}, size(sparse_points,2)); 
                    else
                        color_occupy(free_idx(1)) = 1;
                        fprintf('The to-do-segmentation cluster with label %d and color %s   Traj #: %d \n', i, label_color{free_idx(1)}, size(sparse_points,2)); 
                    end
                else
                    fprintf('The to-do-segmentation cluster with label %d and color %s   Traj #: %d \n', i, label_color{i}, size(sparse_points,2));
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% start processing pipeline
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf('\tComputing motion model\n');
            obj.CurrState = bsmc_inferM(obj.CurrState, obj.Options);
            
            fprintf('\tPredicting Model ');
            tid1=tic;
            % use motion vector to predict app model and label prior 
            % eq(9) p(a | P) and eq(10) p(l = f| P) 
            obj.CurrState = bsmc_predictModelKDE(obj.PrevState, obj.CurrState, obj.Options); % add app and prob model to currstate
            e1 = toc(tid1);
            fprintf('(%f seconds) ...\n', e1);
                
            fprintf('\tComputing Segmentaiton ');
            tid1=tic;
            obj.CurrState = bsmc_computeSegKDE(obj.CurrState, obj.Options); % add debug_info and lmask ,mask % mask is the separation layers of lmask 
            e1 = toc(tid1);
            fprintf('(%f seconds) ...\n', e1);
            
            obj.CurrState.lmaskBfMorph = obj.CurrState.lmask;% lmask will be changed in bsmc_morphProcess. 
            fprintf('\tMorphology process... ');
            obj.CurrState = bsmc_morphProcess(obj.CurrState, obj.Options);
            
            fprintf('\tUpdating Model... ');
            tid1=tic;
            obj.CurrState = bsmc_updateModelKDE(obj.CurrState, obj.Options);
            e1 = toc(tid1);
            fprintf('(%f seconds)\n', e1);  
        end
        
    
        function segIm = genSeg(obj)
            segIm = obj.CurrState.frame;
            if (isfield(obj.CurrState,'seg') && ~isempty(obj.CurrState.seg.mask))
                for j=1:3
                    segIm(:,:,j) = uint8(double(obj.CurrState.frame(:,:,j)) .* obj.CurrState.seg.mask{2} + (double(segIm(:,:,1)) .* ~obj.CurrState.seg.mask{2}) * 0.2);
                end
            end
        end
        
    end
    methods (Static= true)
        function ia = matches_any(set1,set2)
            ia = false(size(set1));
            for i=set2
                ia = ia | (set1 == i);
            end
        end
        
        function im = genAppModel(app, colorspace)
            if strcmp(colorspace, 'rgs')
                im = uint8(bsmc_rgs2RGB(app));
            else
                im = uint8(app);
            end
        end
    end
end
