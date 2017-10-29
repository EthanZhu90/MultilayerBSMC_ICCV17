classdef BroxDataset < handle

    properties
          RootPath = '../Data/moseg_dataset';
    end
    
    methods
        function obj = BroxDataset(rootpath)
            if (nargin == 1)
                obj.RootPath = rootpath;
            else
                % % need to costomize your path to Data folder. 
                obj.RootPath = '../Data/moseg_dataset';
            end            
        end
        
        function genvid(obj)
            seqnames = obj.getSequences();
            for i=1:length(seqnames)
                seq = obj.getSequence(seqnames{i});
                seq.genvid();
            end
        end
        
        function [names, numframes] = getSequences(obj)
            names = dir(obj.RootPath);
            names = names(~cellfun(@(x,y) strcmp(x,'.') | strcmp(x,'..') | ~y,{names(:).name}, {names(:).isdir}));
            names = {names(:).name};
            
            numframes = zeros(1, length(names));
            for i=1:length(names)
                seq = obj.getSequence(names{i});
                numframes(i) = length(seq.Frames);
            end
        end
        
        function seq = getSequence(obj, seqname, numframes)
            if (nargin < 3)
                
                seq = moseg.VideoSeq([obj.RootPath '/' seqname], seqname);
            else
                seq = moseg.VideoSeq([obj.RootPath '/' seqname], seqname, 'numframes', numframes);                
            end
        end
        
        function precomputeOpticalFlow(obj)
            % Precompute optical flow
            seqnames = obj.getSequences();

            for j=1:length(seqnames)
                disp('============================');
                disp([ num2str(j) ': ' seqnames{j}]);
                disp('============================');
                seq = obj.getSequence(seqnames{j});
                seq.precomputeOpticalFlow();
            end
        end
 
    end
end