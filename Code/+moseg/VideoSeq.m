classdef VideoSeq < handle
    properties
        RootPath;
        SeqName;
        SeqRoot;

        Frames;
        ImPrefix;
        ImNumFmt;
        ImExt;
        NumFrames;
        
        % Brox Malik Tracks and Labeling (Lazy Loading)
        Tracks;
        LTracks = [];
        
        % Ground Truth Data (Populated at construction)
        GTFrames;
        GTNumRegions;
        GTScales;
        GTPenalty;
    end
    
    properties (Access = private)
        interlaced = true;
        FileType;
        vr;
    end        
    
    methods
        function obj = VideoSeq(rootpath, seqname, varargin)
            %%% construct function
            obj.RootPath = rootpath;                        
            obj.SeqName = seqname;
            obj.SeqRoot = rootpath;
            
            for i=1:2:length(varargin)
                switch(varargin{i})
                    case 'numframes'
                        numframes = varargin{i+1};
                    case 'interlaced'
                        obj.interlaced = varargin{i+1};
                    otherwise
                        error('Unknown parameter name %s', varargin{i});
                end
            end
                        
            if (exist(rootpath, 'dir'))
                obj.FileType = 'dir';
            else
                obj.FileType = 'file';
                if (exist(rootpath, 'file'))
                    obj.FileType = 'file';
                    obj.vr = cv.VideoCapture(rootpath);
                else
                    error('not a valid file or directory');
                end
            end
            
            if (strcmp(obj.FileType,'dir'))
                % Extract frame numbers and file naming info
                [obj.ImPrefix obj.ImExt] = obj.getImagePrefix();
                [obj.Frames obj.ImNumFmt] = obj.getFrameNumbers();
            else
                if (obj.interlaced)
                    mult = 2;
                else
                    mult = 1;
                end
                obj.Frames = 1:mult*obj.vr.get('FrameCount');
            end
            
            if (exist('numframes','var'))
                obj.Frames = obj.Frames(1:numframes);
            end
            
%             if (strcmp(obj.FileType, 'dir'))
%                 % Read Ground Truth definition file
%                 obj.readDefinition();
%             end
        end
        
        function genvid(obj)
            outpath = '../Videos';
            if (~exist('../Videos', 'dir')); mkdir('../Videos'); end
            vid = VideoWriter([outpath '/' obj.SeqName '.avi']);   
            open(vid);        
            for fno=obj.Frames
                im = obj.readImage(fno);
                writeVideo(vid, im);
            end 
            close(vid);        
        end    
        
        function genvidvr(obj, frames, frmrate)
            if (nargin < 2)
                frames = obj.Frames;
            end
            outpath = '../Videos';
            if (~exist('../Videos', 'dir')); mkdir('../Videos'); end
            vid = VideoWriter([outpath '/' obj.SeqName '-slow.avi']);   
            vid.FrameRate = frmrate;
            open(vid);        
            for fno=frames
                im = obj.readImage(fno);
                writeVideo(vid, im);
            end 
            close(vid);                    
        end
        
        function gentrackingvid(obj,trajlen, frames)
            if (nargin < 2)
                trajlen = 3;
            end
            
            if (nargin < 3)
                frames = obj.Frames;
            end
            
            outpath = '../Videos';
            if (~exist('../Videos', 'dir')); mkdir('../Videos'); end
            vid = VideoWriter([outpath '/' obj.SeqName '-tracking.avi']);   
            open(vid);        
            for fno=frames
                im = obj.gentracking(fno, trajlen);
                writeVideo(vid, im);
                disp(fno)  
            end 
            close(vid);  
        end
        
        function gendetectvid(obj, modelname)
            outpath = '../Videos';
            if (~exist('../Videos', 'dir')); mkdir('../Videos'); end
            vid = VideoWriter([outpath '/' obj.SeqName '-detect.avi']);   
            open(vid);        
            
            pcdir = [obj.SeqRoot '/PrecomputedDetection_' obj.SeqName '/' ];            
            figure(1);
            for j=obj.Frames
                detfile = [pcdir '/' obj.SeqName '_' num2str(j) '_' modelname '_detections.mat'];
                load(detfile);
                im = obj.readImage(j);
            
                imshow(im, 'Border','tight');
                showboxes(im, ds);
                drawnow;                
                f = getframe(gca);
                [im,map] = frame2im(f);    %Return associated image data 
                if isempty(map)            %Truecolor system
                  rgb = im;
                else                       %Indexed system
                  rgb = ind2rgb(im,map);   %Convert image data
                end                            
                writeVideo(vid, rgb);
            end 
            close(vid);  
        end
        
        function rgb = gentracking(obj, imno, trajlen)
            if (nargin < 3)
                trajlen = 3;
            end
            im = obj.readImage(imno);
            
            % Find tracks that have a point at imno and extends to 
            % Find set of tracks which has a point in imno
            StartFrame = obj.Frames(1);
            Tracks = obj.Tracks;
            %Tracks = Tracks(1:10:end);
            TrackStart = cellfun(@(x)x(1), {Tracks(:).frames});
            TrackEnd = cellfun(@(x)x(end), {Tracks(:).frames});
            
            fno = imno - StartFrame + 1;
            tracked = (TrackStart < fno-(trajlen -1) & fno<= TrackEnd);
            tno = find(tracked);
            
            points = zeros(2*trajlen, sum(tracked));
            xx = cellfun(@(x,y) x(1,fno - y(1)-(trajlen-1):fno-y(1)), {Tracks(tracked).points}, {Tracks(tracked).frames},'UniformOutput',false);
            yy = cellfun(@(x,y) x(2,fno - y(1)-(trajlen-1):fno-y(1)), {Tracks(tracked).points}, {Tracks(tracked).frames},'UniformOutput',false);
            points(1:2:end, :) = reshape([xx{:}], trajlen, []);
            points(2:2:end, :) = reshape([yy{:}], trajlen, []);

            figure(1);
            imshow(im, 'Border','tight');
            hold on;
            if (trajlen ==1)
            plot(points(1:2:end,:), points(2:2:end,:), 'g+');
            else
            plot(points(1:2:end,:), points(2:2:end,:), 'g-', 'LineWidth', 3);
            end
            hold off;
            f = getframe(gca);
            [im,map] = frame2im(f);    %Return associated image data 
            if isempty(map)            %Truecolor system
              rgb = im;
            else                       %Indexed system
              rgb = ind2rgb(im,map);   %Convert image data
            end            
        end
              
        function tracks = get.Tracks(obj)
            if (isempty(obj.LTracks))
                obj.LTracks = obj.readtracks();
                tracks = obj.LTracks;
            else
                tracks = obj.LTracks;
            end
        end
        
        function im = readImage(obj, imno)
            if (strcmp(obj.FileType, 'file'))
                if (obj.interlaced)
                    obj.vr.set('PosFrames', floor((imno-1)./2));
                    im = obj.vr.read;
                    if (mod(imno,2) == 0)
                        % Even field
                        im(1:2:end,:,:) = im(2:2:end,:,:);
                    else
                        % Odd field
                        im(2:2:end,:,:) = im(1:2:end,:,:);
                    end
                else
                    obj.vr.set('PosFrames', imno-1);                    
                    im = obj.vr.read;
                end
            else                
                im = imread( sprintf( ...
                    [ obj.SeqRoot '/' obj.ImPrefix obj.ImNumFmt '.' obj.ImExt], imno));
            end
        end        
        
        function im = readGTImage(obj, imno) 
            im = imread(sprintf(...
                [ obj.SeqRoot '/' obj.ImPrefix obj.ImNumFmt '.pgm'],imno));            
        end     
        
        function precomputeOpticalFlow(obj)
            flowpath = [obj.SeqRoot '/OpticalFlow/'];            
            flowmbpath = [obj.SeqRoot '/BroxMalikResults/'];
            for i=obj.Frames(1:end-1)
                disp(i);
                fflowmbfile = [flowmbpath sprintf('ForwardFlow%03d.flo', i-obj.Frames(1))];
                bflowmbfile = [flowmbpath sprintf('BackwardFlow%03d.flo', i-obj.Frames(1))];
                fflowfile = [flowpath sprintf('ForwardFlow%03d.mat', i)];
                bflowfile = [flowpath sprintf('BackwardFlow%03d.mat', i)];

                im1 = obj.readImage(i);
                im2 = obj.readImage(i+1);

                if (~exist(fflowfile, 'file'))
                    if (exist(fflowmbfile,'file'))
                        disp('Found Middleburry file, using.');
                        fflow = moseg.MosegUtils.readMiddlebury(fflowmbfile);                    
                    else                        
                        disp('Cant find file, computing.');                        
                        fflow = moseg.MosegUtils.opticalflow(double(im1), double(im2));
                    end
                    save(fflowfile, 'fflow');
                end
                
                if (~exist(bflowfile, 'file'))
                    if (exist(bflowmbfile, 'file'))
                        disp('Found Middleburry file, using.');                        
                        bflow = moseg.MosegUtils.readMiddlebury(bflowmbfile);                                            
                    else
                        disp('Cant find file, computing.');                                                
                        bflow = moseg.MosegUtils.opticalflow(double(im2), double(im1));
                    end
                    save(bflowfile, 'bflow');
                end
            end
        end
        
        function [tracks, ncluster] = readtracks(obj)
            trackfile = dir(sprintf('%s/BroxMalikResults/Tracks%d.dat',obj.SeqRoot, length(obj.Frames)));
            fname = sprintf('%s/BroxMalikResults/%s',obj.SeqRoot,trackfile.name);
            fid = fopen(fname,'r');
            aLength = fscanf(fid,'%d',1);
            aTrackNo = fscanf(fid,'%d',1);
            
            % Create a structure Array
            points = cell(1,aTrackNo);
            frames = cell(1,aTrackNo);
            labels = cell(1,aTrackNo);
            
            for i = 1:aTrackNo 
                labels{i} = fscanf(fid,'%d',1) + 1;
                aSize = fscanf(fid,'%d',1);
                points{i} = zeros(2,aSize);
                frames{i} = zeros(1,aSize);
                data = fscanf(fid,'%f %f %d', [3 aSize]);
                points{i}= data(1:2,:)+1; % Convert to 1-indexed coordinates
                frames{i} = data(3,:)+1;
            end
            ncluster = 0;
            for i=1:aTrackNo
                if labels{i} > ncluster
                    ncluster = labels{i}+1;
                end
            end
            fclose(fid);
            
            tracks = struct('points', points, 'frames', frames, 'label', labels);
        end
        
        function fflow = getForwardFlow(obj, imno)
            flowmbpath = [obj.SeqRoot '/BroxMalikResults/'];
            fflowmbfile = [flowmbpath sprintf('ForwardFlow%03d.flo', imno-obj.Frames(1))];

            if (exist(fflowmbfile,'file'))
                disp('Found Middleburry file, using.');
                fflow = moseg.MosegUtils.readMiddlebury(fflowmbfile);                    
            else                        
                error('Cant find forward flow file for image %d of sequence %s.', imno, obj.SeqName);                        
            end
        end    
        
        function bflow = getBackwardFlow(obj, imno)
            flowmbpath = [obj.SeqRoot '/BroxMalikResults/'];
            bflowmbfile = [flowmbpath sprintf('BackwardFlow%03d.flo', imno-obj.Frames(1))];

            if (exist(bflowmbfile, 'file'))
                disp('Found Middleburry file, using.');                        
                bflow = moseg.MosegUtils.readMiddlebury(bflowmbfile);                                            
            else
                error('Cant find file, computing.');                                                
            end
        end
    end
    
    methods (Access = private)
        

        
        function [tmpl, suff] = getImagePrefix(obj)            
            % Find the number of frame automaticaly
            files = dir([obj.SeqRoot '/*.ppm']);
            if (isempty(files))
                suff = 'jpg';  %suff = 'png'; by Eason
                files = dir([obj.SeqRoot '/*.jpg']);  %files = dir([obj.SeqRoot '/*.png']); by Eason            
            else
                suff = 'ppm';
            end
            m = regexp(files(1).name,['([-a-zA-Z0-9]+_*)(\d*).' suff],'tokens');                
            if (isempty(m{1}{2}))
                m = regexp(files(1).name,['([-a-zA-Z]+_*)(\d*).' suff],'tokens');                
            end
            tmpl = m{1}{1}; 
        end
                
        function [frames numformat] = getFrameNumbers(obj)
            % Find the number of frame automaticaly
            files = dir([obj.SeqRoot '/*.' obj.ImExt]);
            last_frame = 0;
            minlen = 1;
            maxlen = 100;
            start_frame = 1000000;
            for i=1:length(files)
                m = regexp(files(i).name,['(' obj.ImPrefix '_*)(\d*).' obj.ImExt],'tokens');                
                minlen = max(minlen, length(m{1}{2}));
                maxlen = min(maxlen, length(m{1}{2}));
                frameNo = str2num(m{1}{2}); %sscanf(files(i).name, [seqname '%03d.ppm']);
                last_frame = max(last_frame, frameNo);
                start_frame = min(start_frame, frameNo);
            end            
            frames = start_frame:last_frame;
            
            if (minlen == maxlen)
                numformat = ['%0' num2str(minlen) 'd'];
            else
                numformat = '%d';
            end
            
        end        
    end
end