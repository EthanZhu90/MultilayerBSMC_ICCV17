classdef Tracker < handle
    methods 
        function SetSequence(obj, seq)
        end
    end
    methods(Abstract)
        [points, tracked] = step(obj, im1, im2, imno, points);        
    end
end