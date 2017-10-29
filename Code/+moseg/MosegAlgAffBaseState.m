classdef MosegAlgAffBaseState < handle
    properties 
        memTrackno = [];
        memTrajIds = [];      % Set of tracks in embedding space.                
        maxTrajId = 0;
        memStartFrame = [];
        memEndFrame = [];
        memLbls = [];
        memProb = [];
        frameNo;
        lowAff = [];
        outliers = [];
        bgclust = 0;
    end
    
    properties(Transient)
        W = [];
        Mask = [];
        Dm = [];
        Dsp = [];     
        DspInRange =[]; 
    end
    
    methods
        function obj =MosegAlgAffBaseState(fno)
            obj.frameNo = fno;            
        end
    end
end
