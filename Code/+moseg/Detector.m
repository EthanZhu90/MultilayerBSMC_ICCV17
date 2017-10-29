classdef Detector < handle
    methods 
        function SetSequence(obj, seq)
        end
    end
    methods(Abstract)
        [points, trackno] = detect(obj, ~, p_in, imno);
    end    
end