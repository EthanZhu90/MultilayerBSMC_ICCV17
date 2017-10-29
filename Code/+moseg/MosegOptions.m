classdef MosegOptions
    properties
        lookahead = 5;
        buffSize = 5; % Number of points to hold in memory for each trajectory.
        dropEnded = true;

        % Options for tracking.
        maxTrajFrame = 50000; 
        
        % Maintianing extra datastructures (slows things down).
        extrads = false;
    end
    
    methods
        function obj = MosegOptions(buffSize)
            if (nargin > 0)
                obj.buffSize = buffSize;
            end
        end
    end
end