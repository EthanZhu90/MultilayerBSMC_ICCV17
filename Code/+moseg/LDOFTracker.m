classdef LDOFTracker < moseg.Tracker
    properties 
        FlowPath
    end
    methods        
        
        function obj = LDOFTracker(flowPath)            
            if (nargin > 0)
                if (~exist(flowPath, 'dir'))
                    mkdir(flowPath);
                end

                obj.FlowPath = flowPath;
            end
        end
        
        function SetSequence(obj, seq)
            obj.FlowPath = [seq.SeqRoot '/OpticalFlow/'];  
            if (~exist(obj.FlowPath, 'dir'))
                    mkdir(obj.FlowPath);
            end
            
        end
        
        function [fflow, bflow] = getFlow(obj, im1, im2, imno)
            % Compute optical flows between frames imno-1 and i
            fflowfile = [obj.FlowPath sprintf('ForwardFlow%03d.mat', imno)];
            %disp(fflowfile);
            if (exist(fflowfile, 'file'))
                %disp('Loading precomputed Forward Optical flow');                    
                load(fflowfile, 'fflow');
            else
                %disp('Computing Forward Optical flow');            
                fflow = moseg.MosegUtils.opticalflow(double(im1), double(im2));
                save(fflowfile, 'fflow');
            end

            bflowfile = [obj.FlowPath sprintf('BackwardFlow%03d.mat', imno)];
            if (exist(bflowfile, 'file'))
                %disp('Loading precomputed Backward Optical flow');            
                load(bflowfile, 'bflow');
            else
                %disp('Computing Backward Optical flow');    
                bflow = moseg.MosegUtils.opticalflow(double(im2), double(im1));
                save(bflowfile, 'bflow');
            end            
        end
        
        function [points, tracked] = step(obj, im1, im2, imno, points)
            [fflow, bflow] = obj.getFlow(im1, im2, imno);
            [points, tracked] = oftrack(fflow, bflow, points);    
        end
    end
end