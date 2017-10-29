% Just assign all labels the same class.
classdef MosegAlgDummy < moseg.MosegAlg
    properties 
        Name = 'Dummy';
        Options = struct('dataset', 'seqname');
    end
    methods
        function obj = MosegAlgDummy()            
            obj.SupportsPartialTracks = true;
        end
    end
    
    methods (Access = protected)
        function lbls = stepPartial(obj, ~, ~, trackno, ~, ~)
            lbls = ones(1,length(trackno));
        end
    end
end