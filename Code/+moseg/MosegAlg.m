classdef MosegAlg < handle
    properties (Abstract = true)
        Name;        
        Options;        
    end

    properties
        SupportsPartialTracks = false;
        TimingStats = struct;
    end
    
    
    methods
        function lbls = step(obj, points, trajIds, trackno, len, frameNo)
            if (~obj.SupportsPartialTracks)
                lbls = zeros(1,size(points,2));                        
                fulltraj = len >= obj.Options.windowSize;
                if (sum(fulltraj) ~=0)
                    lbls(fulltraj) = obj.stepFull(points(1:2*obj.Options.windowSize,fulltraj), ...
                        trajIds(fulltraj), ...
                        trackno(fulltraj), ...
                        len(fulltraj), ...
                        frameNo);
                end
            else
                lbls = obj.stepPartial(points, ...
                    trajIds, ...
                    trackno, ...
                    len, ...
                    frameNo);                
            end
        end
        
        function SetSequence(obj, seq)
            % Function for Algorithms that needs the seqeunce
        end
        
        function visualize(obj)
        end
        
        function printStats(obj)
            % Do nothing here.
            fprintf('\n');
        end
        
        function reset(obj, startframe)
            % If the algorithm is statefull then reset datastructures.
        end
    end
    
    methods (Access = protected)
        function lbls = stepFull(obj, points, trajIds, trackno, len, frameNo)
            lbls = zeros(1,size(points,2));
        end
        function lbls = stepPartial(obj, points, trajIds, trackno, len, frameNo)
            lbls = zeros(1,size(points,2));            
        end
        function debugmsg(obj, str, varargin)
            if (obj.Options.verbose)
                fprintf(str, varargin{:});
            end
        end        
        
    end
    
end