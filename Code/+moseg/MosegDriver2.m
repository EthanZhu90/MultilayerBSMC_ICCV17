classdef MosegDriver2 < handle    
    properties 
        mosegmenter;
    end
    
    properties (Access = public)  %private changed by Ethan 10/23
        curframe;
      
        % Tracking stuff
        points;
        trackno;
        maxTrackno;
        tracker;
        detector;
        im1;

        % Motionsegmentation
        alg;
        startframe;
    end
    
    methods
        function obj = MosegDriver2(varargin)
            
            if (mod(length(varargin),2) ~= 0)
                error('Incorrect number of arguments');
            end
            
            % Setup defaults
            obj.startframe = 1;
            obj.detector = moseg.GridDetector;
            %obj.tracker = moseg.BidirectionalTracker(moseg.LKTracker);
            obj.tracker = moseg.LDOFTracker();
            obj.alg = moseg.MosegAlgLblProp(obj.startframe);
            
            % Hold at most 5 points            
            buffSize = 5;
            for i=1:2:length(varargin)
                switch(varargin{i})
                    case 'Algorithm'
                        obj.alg = varargin{i+1};
                    case 'Detector'
                        obj.detector = varargin{i+1};
                    case 'Tracker'
                        obj.tracker = varargin{i+1};
                    case 'StartFrame'
                        obj.startframe = varargin{i+1};
                    case 'BuffSize'
                        buffSize = varargin{i+1};
                    otherwise
                        error('Unrecognized parameter name %s', varargin{i});
                end
            end
                        
            % Reset the algorithm
            obj.alg.reset(obj.startframe); % For stateful algorithms.            
            obj.curframe = obj.startframe-1;
            
            opt = moseg.MosegOptions(buffSize); 
            obj.mosegmenter = moseg.Moseg(obj.alg, opt);
            
            % Initialize data structures
            obj.points = zeros(2,0);
            obj.trackno = [];
            obj.maxTrackno = 0;                     
            
        end
        
        function step(obj, im2)
            import moseg.*;       
            tid= tic;
            obj.curframe = obj.curframe + 1;
            
            fprintf('============================\n');
            fprintf('Processing Frame %03d (%d).\n', obj.curframe, obj.curframe-obj.startframe+1);
            
            
            tid2 = tic;
            
            ntracked = 0;
            if (obj.curframe ~= obj.startframe)
                % Track points
                [obj.points, tracked] = obj.tracker.step(obj.im1, im2, obj.curframe-1, obj.points);
                obj.points = obj.points(:, tracked);
                obj.trackno = obj.trackno(tracked);
                ntracked = length(obj.points);
            end
            tracktime = toc(tid2);
            
            % Detect new points
            newpoints = obj.detector.detect(im2, obj.points, obj.curframe);
            newtrackno = obj.maxTrackno + (1:size(newpoints,2));
            obj.maxTrackno = obj.maxTrackno + size(newpoints,2);
            
            % Append tracked points to detected points
            obj.trackno = [obj.trackno newtrackno];
            obj.points = [obj.points newpoints];

            %obj.alg.Frame = im2;
            obj.mosegmenter.step(obj.points, obj.trackno);            
       
            obj.im1 = im2;
            
            elapsedtime = toc(tid);
            fprintf('\tCompleted Frame %03d (%d) . PointsTracked %d. Tracks: %d. Elapsed time %f\n', ...
                obj.curframe, obj.curframe-obj.startframe+1, ntracked, size(obj.points,2), elapsedtime);

            fprintf('\tTiming Stats: tracktime=%f, ', tracktime); 
            obj.mosegmenter.printStats();
            fprintf('\tAlgorithm Stats: ');        
            obj.alg.printStats();
        end

        function lookahead(obj, im2)
            import moseg.*;       
            tid= tic;
            obj.curframe = obj.curframe + 1;
            
            fprintf('============================\n');
            fprintf('Processing Frame %03d (%d).\n', obj.curframe, obj.curframe-obj.startframe+1);
            
            
            tid2 = tic;
            
            ntracked = 0;
            if (obj.curframe ~= obj.startframe)
                % Track points
                [obj.points, tracked] = obj.tracker.step(obj.im1, im2, obj.curframe-1, obj.points);
%                 [obj.points, tracked, ttrackno] = obj.tracker.step(obj.im1, im2, obj.curframe-1, obj.points);
                obj.points = obj.points(:, tracked);
                obj.trackno = obj.trackno(tracked);
%                 assert(all(ttrackno(tracked) == obj.trackno));
                ntracked = length(obj.points);
            end
            tracktime = toc(tid2);
            
            % Detect new points
            newpoints = obj.detector.detect(im2, obj.points, obj.curframe);
            newtrackno = obj.maxTrackno + (1:size(newpoints,2));
            obj.maxTrackno = obj.maxTrackno + size(newpoints,2);
            
            % Append tracked points to detected points
            obj.trackno = [obj.trackno newtrackno];
            obj.points = [obj.points newpoints];

            %obj.alg.Frame = im2;
            obj.mosegmenter.lookahead(obj.points, obj.trackno);            
            
            obj.im1 = im2;
            
            elapsedtime = toc(tid);
            fprintf('\tCompleted Frame %03d (%d) . PointsTracked %d. Tracks: %d. Elapsed time %f\n', ...
                obj.curframe, obj.curframe-obj.startframe+1, ntracked, size(obj.points,2), elapsedtime);

            fprintf('\tTiming Stats: tracktime=%f, ', tracktime); 
            obj.mosegmenter.printStats();
            fprintf('\tAlgorithm Stats: ');        
            obj.alg.printStats();
        end
        
        % Simple Ploting, Should useing MosegGUI for interactivity
        function plot(obj, im, varargin)
            ha = gca;
            lookahead = 0;
            for i=1:2:length(varargin)
                switch(varargin{i})
                    case 'Parent'
                        ha = varargin{i+1};
                    case 'LookAhead'
                        lookahead = varargin{i+1};                        
                    otherwise
                        error('Unknown parameter name %s', varargin{i});
                end
            end
            
            im = obj.mosegmenter.genim(im, true, 'LookAhead', lookahead);
            imshow(im, 'Border','tight','InitialMagnification','fit', 'Parent',ha );            
            
            drawnow;            
        end        

    end
    
    methods (Static=true)
           
        function printstats(stats)
            fprintf('Evaluation results for %s\n', stats.seqname);
            fprintf('MoSegEval Version 1.0\n');
            fprintf('Number of frames used from the sequences:\n');
            fprintf('%d\n', stats.numframes);
            fprintf('Number of labeled frames in this time window:\n');
            fprintf('%d\n', stats.numgtframes);
            fprintf('--------------------------\n');
            fprintf('Density (in percent):\n');
            fprintf('%f\n', stats.density);
            fprintf('--------------------------\n');
            fprintf('Overall (per pixel) clustering error (in percent):\n');
            fprintf('%f\n', stats.overallerr);
            fprintf('Clustering error per region (in percent):\n');
            fprintf('Region %d: \n%f\n',[0:stats.nregions-1; stats.rgnerr]);
            fprintf('Visible regions in the evaluated part of the shot:\n');
            fprintf('%d\n', stats.visrgns);
            fprintf('--------------------------\n');
            fprintf('Average (per region) clustering error (in percent):\n');
            fprintf('%f\n', stats.avgerr);
            fprintf('--------------------------\n');
            fprintf('Number of clusters merged to obtain this result (oversegmentation error):\n');
            fprintf('%d\n', stats.oversegpenalty);
            fprintf('--------------------------\n');
            fprintf('Number of regions with less than 10%% error (excluding background):\n');
            fprintf('%d\n', stats.lt10percent);            
        end
        
    end
end