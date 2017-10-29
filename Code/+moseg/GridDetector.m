classdef GridDetector < moseg.Detector
    properties
        stepsize;
    end
    
    methods
        function obj = GridDetector(stepsize)
            if (nargin <1)
                stepsize = 4;
            end
            obj.stepsize = stepsize;
        end
        
        function points = detect(obj, im, points, ~)
            if (nargin < 3);
                points = zeros(2,0);
            end
            if (exist('griddetect','file') == 0)
                addpath('../Utils/OFTrack'); 
            end
            points = griddetect(im, obj.stepsize, points);
        end
    end
end
