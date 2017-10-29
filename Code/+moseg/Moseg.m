classdef Moseg < handle
    % This class follows the strategy design pattern.
    properties
        % Model
        Options = moseg.MosegOptions;
        CurrState;
        MaxLbl = 0;
        
        % View
        Highlight;
        
        TimingStats;
    end
    
    properties(Transient=true)
        PrevState;        
        MosegAlg;
    end
    
    
    methods
        function obj = Moseg(mosegAlg, opt)
            if (nargin > 1)
                obj.Options = opt;
            end

            obj.CurrState.len = zeros(1,0);
            obj.CurrState.points = nan(2*obj.Options.buffSize,0);
            obj.CurrState.trajIds = [];
            obj.CurrState.trackno = [];
            obj.CurrState.lbls = [];
            obj.CurrState.maxTrajId = 0;
            obj.CurrState.startFrame = [];
            obj.CurrState.endFrame = [];
            obj.CurrState.frameNo = 0;
            obj.CurrState.compTime = 0;
            obj.CurrState.totalTime = 0;
            obj.CurrState.bgclust = [];
            obj.CurrState.memTrajIds = []; %%% add by ethan 10/13

            obj.MosegAlg = mosegAlg;
            
        end

        
        function step(obj, points, trackno)
            obj.stepInit();
            obj.stepTraj(points, trackno);
            assert(size(obj.CurrState.points,2) == length(obj.CurrState.trajIds));
            
            tid = tic;
            % Use the set of points that span the entire window.
            % the main process, the main reason of time-consuming 
            obj.CurrState.lbls = obj.MosegAlg.step(obj.CurrState.points, ...
                obj.CurrState.trajIds, ...
                obj.CurrState.trackno, ...
                obj.CurrState.len, ...
                obj.CurrState.frameNo);
            
            M = mode(obj.CurrState.lbls(:));

            if M ~= 0 && isempty(obj.CurrState.bgclust)
                obj.CurrState.bgclust = M;
                fprintf('\n##### Decide the bglabel is %d\n',  obj.CurrState.bgclust);
            elseif( M ~= 0 && ~isempty(obj.CurrState.bgclust))
                bg_ratio = sum(obj.CurrState.lbls(:) == obj.CurrState.bgclust)/numel(obj.CurrState.lbls);
                if bg_ratio < 0.05
                    fprintf('\n##### Switch the bglabel from %d to %d\n',  obj.CurrState.bgclust, M);
                    obj.CurrState.bgclust = M;
                end
            end
            obj.CurrState.compTime = toc(tid);
            obj.CurrState.totalTime = obj.CurrState.totalTime + obj.CurrState.compTime;
        end
        

        
        function printStats(obj)
            s = fieldnames(obj.MosegAlg.TimingStats);
            if (~isempty(s))
                for i=1:length(s)-1
                    fprintf('%s: %f, ', s{i}, obj.MosegAlg.TimingStats.(s{i}));
                end
                fprintf('%s: %f\n', s{end}, obj.MosegAlg.TimingStats.(s{end}));                
            end
        end

        
    end
    
    methods (Access = private)
        function stepInit(obj)
            obj.PrevState = obj.CurrState;
            obj.CurrState = struct;
            obj.CurrState.len = zeros(1,0);
            obj.CurrState.points = nan(2*obj.Options.buffSize,0);
            obj.CurrState.trajIds = [];
            obj.CurrState.lbls = [];
            obj.CurrState.trackno = [];
            obj.CurrState.maxTrajId = 0;
            obj.CurrState.frameNo = obj.PrevState.frameNo + 1;
            obj.CurrState.compTime = 0;
            obj.CurrState.totalTime = obj.PrevState.totalTime;
            obj.CurrState.bgclust = obj.PrevState.bgclust;
            %obj.CurrState.frame = frame;

            if (obj.Options.extrads)
                obj.CurrState.startFrame = [];
                obj.CurrState.endFrame = [];
            end
        end
        
        function stepTraj(obj, dpoints, trackno)
            
            %%% Match points to previous frame
            if (obj.PrevState.maxTrajId > 0)
                [~,ia,ib] = intersect(trackno, obj.PrevState.trackno);
                tracked = false(size(obj.PrevState.points,2),1);
                tracked(ib) = 1;
                points = nan(2,size(obj.PrevState.points,2));
                points(:,ib) = dpoints(:,ia);

                points = [ points ; obj.PrevState.points(1:(obj.Options.buffSize*2-2),:)];

                % Copy and update prevState trajectories
                obj.CurrState.len = obj.PrevState.len; 
                obj.CurrState.len(tracked) = obj.CurrState.len(tracked) + 1;
                obj.CurrState.points = points;
                obj.CurrState.lbls = obj.PrevState.lbls;
                obj.CurrState.trajIds = obj.PrevState.trajIds;
                obj.CurrState.trackno = obj.PrevState.trackno;

                if (obj.Options.extrads)
                    obj.CurrState.startFrame = obj.PrevState.startFrame;
                    obj.CurrState.endFrame   = obj.PrevState.endFrame;
                    obj.CurrState.endFrame(tracked) = obj.CurrState.frameNo;
                end                
                
                if (obj.Options.dropEnded)
                    obj.CurrState.len = obj.CurrState.len(tracked);
                    obj.CurrState.points = obj.CurrState.points(:,tracked);
                    obj.CurrState.lbls = obj.CurrState.lbls(tracked);
                    obj.CurrState.trajIds = obj.CurrState.trajIds(tracked);
                    obj.CurrState.trackno = obj.CurrState.trackno(tracked);
                    
                    if (obj.Options.extrads)
                        obj.CurrState.startFrame = obj.CurrState.startFrame(tracked);
                        obj.CurrState.endFrame   = obj.CurrState.endFrame(tracked);
                    end                                    
                end
            end

            % Detect new points in current frame and add them to set of points
            [~,ia] = setdiff(trackno, obj.PrevState.trackno);
            num_select = length(ia);
            new_points = [dpoints(:,ia); nan(2*(obj.Options.buffSize-1),num_select)];
            new_trackno = trackno(ia);
            
            % Update Datastructures
            obj.CurrState.len = [obj.CurrState.len ones(1,num_select)];
            obj.CurrState.points = [obj.CurrState.points new_points];
            obj.CurrState.trackno = [obj.CurrState.trackno new_trackno];            
            obj.CurrState.trajIds = [obj.CurrState.trajIds (obj.PrevState.maxTrajId+1:obj.PrevState.maxTrajId+size(new_points,2))];
            obj.CurrState.maxTrajId = obj.PrevState.maxTrajId + size(new_points,2);
            obj.CurrState.lbls = [obj.CurrState.lbls zeros(1, num_select)];

   
            % Version that holds start and ending frames for all trajectories.
            if (obj.Options.extrads)
                obj.CurrState.startFrame = [obj.CurrState.startFrame repmat(obj.CurrState.frameNo,1, size(new_points,2))];
                obj.CurrState.endFrame   = [obj.CurrState.endFrame repmat(obj.CurrState.frameNo,1, size(new_points,2))];
                %assert(length(obj.CurrState.startFrame) == obj.CurrState.maxTrajId && length(obj.CurrState.endFrame) == obj.CurrState.maxTrajId);
            end
        end    
    end
    

end

