classdef MosegUtils
    methods(Static)
        function writetracks(tracks, fname) 
            %trackfile = dir(sprintf('%s/BroxMalikResults/Tracks%d.dat',obj.SeqRoot, length(obj.Frames)));
            %fname = sprintf('%s/BroxMalikResults/%s',obj.SeqRoot,trackfile.name);
            fid = fopen(fname,'w');
            % Find last frame
            aLength = max([tracks{:}.frames]);
            fprintf(fid,'%d\n',aLength);
            fprintf(fid,'%d\n',length(tracks));
                        
            for i = 1:length(tracks)
                fprintf(fid,'%d ',tracks{i}.label);
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
        
        function im = latex2png(str)
            prefix = '\documentclass[preview]{standalone} \begin{document}';
            suffix = '\end{document}';
            str = [prefix str suffix];
            fid = fopen('/tmp/im.tex','w+');
            fprintf(fid, '%s', str);
            fclose(fid);
            system('pdflatex /tmp/im.tex');
            system('convert -density 300 im.pdf -quality 90 im.png');
            im = imread('im.png');
        end
        
        function str = latextable(X,varargin)
            % LATEXTABLE produces a table for usage with LaTeX from a numeric array
            % Copyright (c) 2009, Andrew E. Slaughter
            % All rights reserved.
            %
            % Redistribution and use in source and binary forms, with or without
            % modification, are permitted provided that the following conditions are met:
            %     * Redistributions of source code must retain the above copyright
            %       notice, this list of conditions and the following disclaimer.
            %     * Redistributions in binary form must reproduce the above copyright
            %       notice, this list of conditions and the following disclaimer in the
            %       documentation and/or other materials provided with the distribution.
            %     * Neither the name of the <organization> nor the
            %       names of its contributors may be used to endorse or promote products
            %       derived from this software without specific prior written permission.
            %
            % THIS SOFTWARE IS PROVIDED BY <copyright holder> ''AS IS'' AND ANY
            % EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
            % WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
            % DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
            % DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
            % (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
            % LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
            % ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
            % (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
            % SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
            %__________________________________________________________________________
            % SYNTAX:
            %   latextable(X);
            %   latextable(X,'PropertyName',PropertyValue,...)
            %
            % DESCRIPTION:
            %   latextable(X) produces a gerneic table without headers and asks the
            %       for a name and location for saving the file, X must be an 2-D
            %       array of numbers or a cell array, see below.
            %   latextable(X,'PropertyName',PropertyValue,...) allows for custom
            %       inputs to be used for producing the table.
            %
            % INPUT "X":
            % The input may be a purely numeric array or a cell array of mixed
            % components.  The following examples illustrate the type of input
            % accepted.
            %   >> X = [1,2; 3,4];
            %   >> X = {1,'2';'One',4};
            %
            % FUNCTION PROPERTIES:
            %   'Horiz' - Cell array containing labels for horizontal direction, the
            %       cell array must have the same horizontal dimension (colums) as the
            %       "X" input but can contain one or many rows. The default is an empty
            %       cell, thus no labels.
            %   'Vert' - Cell array containing labels for vertical direction, the
            %       cell array must have the same vertical dimension (rows) as the
            %       "X" input but can contain one or many columns. The default is an
            %       empty cell, thus no labels.
            %   'Vline' - numeric array prescribing locations of vertical lines, the
            %       line is inserted after the specified column. For example, using
            %       [1,2] inserts lines after the 1st and 2nd columns. A value of NaN
            %       at the end ([...,NaN]) places a line at the end of the table. The
            %       default is [], no vertical lines.  The column numbers include the
            %       labels as columns.
            %   'Hline' - numeric array prescribing locations of horizontal lines,
            %       which operates similar to 'Vline', the default is [0,NaN]
            %   'format' - specifies the numeric format, e.g. %3.1f.
            %   'name' - string specifing the filename
            %   'save' - string specifing the name to save the supplied setup
            %   'load' - the setup to be loaded, all properties prior to usage are
            %       overwritten, those after are applied.  Note, the 'default' is
            %       automatically used.  To view the available settings, type the
            %       following into the command-line: "getpref('latextable_settings')"
            %
            % EXAMPLES:
            %   >> latextable(rand(3,3),'Horiz',{'1','2','3'},'Hline',[0,1,NaN]);
            %   >> latextable(rand(3,3),'Hline',[],'save','nolines');
            %   >> latextable(rand(5,5),'load','nolines');
            %
            % PROGRAM OUTLINE:
            % 1 - INITILIZE PROGRAM
            % 2 - CONVERT INPUT INTO A SINGLE COLUMN OF '&' DELIMINATED STRINGS
            % 3 - ADD VERTICAL HEADINGS TO STRINGS
            % 4 - BUILD HEADER ROW(S)
            % 5 - COMBINE ROWS AND ADD HORIZONTAL LINES
            % 6 - BUILD COLUMN INSTRUCTIONS
            % 7 - OUTPUT DATA TO FILE
            % 8 - APPLY PREFERENCES
            % SUBFUNCTION: subfunction_options
            % SUBFUNCTION: convertdata
            % APPPYPREF applies the last used directory and property settings
            %__________________________________________________________________________
            
            % 1 - INITILIZE PROGRAM
            a = subfunction_options(varargin);
            
            % 2 - CONVERT INPUT INTO A SINGLE COLUMN OF '&' DELIMINATED STRINGS
            r = size(X,1); c = size(X,2);
            A = convertdata(X,a.format);
            
            % 3 - ADD VERTICAL HEADINGS TO STRINGS
            Hprfx = '';
            if ~isempty(a.vert) && size(a.vert,1) == r;
                for i = 1:r
                    A{i,1} = [sprintf('%s & ',a.vert{i,:}),A{i}];
                end
                Hprfx(1:size(a.vert,2)) = '&';
            end
            
            % 4 - BUILD HEADER ROW(S)
            H = {};
            if ~isempty(a.horiz)&& size(a.horiz,2) == c ;
                for i = 1:size(a.horiz,1)
                    H{i,1} = sprintf('%s & ',a.horiz{i,:});
                    H{i} = [Hprfx,' ',regexprep(H{i},'& ','\\\',c)];
                end
            end
            
            % 5 - COMBINE ROWS AND ADD HORIZONTAL LINES
            ROWS = [H;A]; h = a.hline; nr = length(ROWS);
            for i = 1:length(h);
                if h(i) == 0;
                    ROWS = ['\hline'; ROWS(1:nr)];
                    h = h + 1; nr = nr + 1;
                elseif isnan(h(i));
                    ROWS{nr+1} = '\hline';
                    nr = nr + 1;
                else
                    ROWS = [ROWS(1:h(i));'\hline';ROWS(h(i)+1:nr)];
                    h = h + 1; nr = nr + 1;
                end
            end
            
            % 6 - BUILD COLUMN INSTRUCTIONS
            col(1:c+size(a.vert,2)) = 'c';
            v = a.vline;
            for i = 1:length(v);
                if v(i) == 0; col = ['|',col];
                elseif isnan(v(i)); col = [col,'|'];
                else
                    col = regexprep(col,'c','c|',v(i));
                end
            end
            
            % 7 - OUTPUT DATA TO FILE
            str = sprintf('%s\n',['\begin{tabular}{',col,'}']);
            for i = 1:length(ROWS);
                str = [str sprintf('%s\n',ROWS{i})];
            end
            str = [str sprintf('%s\n','\end{tabular}')];
            
            %--------------------------------------------------------------------------
            % SUBFUNCTION: subfunction_options
            function a = subfunction_options(in)
                % SUBFUNCTION_OPTIONS seperates user inputed property modifiers
                
                % 1 - SET THE DEFAULT SETTINGS
                % 1.1 - Standard settings
                a.horiz  = {};          a.vert  = {};
                a.hline  = [0,NaN];     a.vline = [];
                a.format = '%3.1f';
                
                % 2 - SEPERATE THE DATA FROM OPTIONS
                list = fieldnames(a); n = length(in);
                for k=1:2:n
                    opt = in{k}; value = in{k+1}; match = '';
                    
                    % 2.2 - Compare modifier tag with available options
                    if ischar(opt);
                        match = strcmpi(opt,list);
                        if ~isempty(match);
                            a.(list{match}) = value;
                        end
                    end
                    
                    % 2.4 - Produce an error message if modifier is not found
                    if isempty(match);
                        disp( ['The property modifier, ',opt,', was not recoignized.']);
                    end
                end
                
                % 4 - CORRECT VERTICAL AND HORIZONTAL ARRAYS WITH A SINGLETON DIMENSION
                if size(a.horiz,2) == 1; a.horiz = a.horiz'; end
                if size(a.vert,1) == 1; a.vert = a.vert'; end
            end
            
            function A = convertdata(X,form)
                % CONVERTDATA changes numeric and cell arrays into strings for each row
                % X must be a numeric array or a cell array, with a string of numeric value
                % in each entry.
                
                % Define the size of the array
                r = size(X,1); c = size(X,2);
                
                % For completely numeric arrays
                if isnumeric(X);
                    for j = 1:r;
                        A{j,1} = sprintf([form,' & '],X(j,:));
                        A{j} = regexprep(A{j},'& ','\\\',c);
                    end
                    
                    % For cell arrays
                elseif iscell(X);
                    for k = 1:r;
                        A{k,1} = num2str(X{k,1},form);
                        for j = 2:c;
                            A{k,1} = [A{k,1},' &',num2str(X{k,j},form)];
                        end
                        A{k,1} = [A{k,1},'\\'];
                    end
                end
            end
            
        end
        
        function str = exporttxt(X)
            str = [];
            line = [];
            varnames = fieldnames(X)';
            line2= sprintf([repmat('=', 1, 20*length(varnames)) '\n']);
            for j=varnames
                line = [line sprintf('%-20s', j{1})];                
            end
            line = [line sprintf('\n')];
            str = [str line line2];

            for i=1:length(X)
                line = [];
                for j=varnames
                    val = X(i).(j{1});                        
                    if (ischar(X(i).(j{1})))
                        fmt = '%-20s';
                    else
                        fmt = '%-20g';
                        if (length(val) ~=1)
                            val = sprintf('%g,', val);
                            fmt = '%-20s';
                        end
                    end
                    line = [line sprintf(fmt, val)];
                end
                line = [line sprintf('\n')];
                str = [str line];
            end
            str = [str line2];
        end
        
        function colors = gencolormap(scale)
            if (nargin < 1)
                scale = true;
            end
            
            colorsa = [ 1 0 0;
                0.75 0.25 0;
                0.5 0.5 0;
                0.25 0.75 0;
                0 1 0];

%             colorsa = [ 1 0 0;
%                 0 1 0];

            colors = ...
                [colorsa;
                colorsa(:,[3 1 2]);
                colorsa(:,[2 3 1])];
            
            if (scale)
                colors = round(255 * colors);
            end
            
            % Reshuffle the colors
            idx = reshape(1:(size(colorsa,1)*3), size(colorsa,1),3);
            idx = idx';
            idx = idx(:);
            colors = colors(idx,:);
            colors = repmat(colors, 20,1);
        end
        
        function flow = opticalflow(im1, im2)
            if (exist('mex_LDOF','file') == 0)
                addpath('../Utils/pami2010Matlab');
            end
            flow = mex_LDOF(double(im1),double(im2));
        end
        
        function [trackedx, trackedno, detectx] = detectandtrack(im1, im2, flow1, flow2, x, trackno)
            if (exist('ldoftrack','file') == 0)
                addpath('/home/elqursh/Projects/Research/Libraries/Tracking/trackingLinux64');
            end
            if (nargin < 5)
                x = zeros(2,0);
            end
            disp('Detecting features and tracking');
            [xout, lbl] = ldoftrack(im1, im2, flow1, flow2, x);
            n = size(x,2);
            trackedlbl = lbl(1:n);
            detectedlbl = lbl((n+1):end);
            
            trackedx = xout(:,1:n);
            trackedx = trackedx(:, trackedlbl);
            trackedno = trackno(trackedlbl);
            
            detectx = xout(:,(n+1):end);
            detectx = detectx(:, detectedlbl);
        end
        
        function Wnn = nnSparseMat(W, K_nn, K_random, thresh, binary)
            [sorted,index] = sort(W,'descend');
            n=length(W);
            
            % Since we know the number of nearest neighbors we can
            % preallocate matrices
            ix_nn = repmat(1:n,K_nn,1);
            jx_nn = index(2:(1+K_nn),:);
            if (binary)
                v_nn = ones(K_nn,n);
            else
                v_nn = sorted(2:(1+K_nn),:);
            end
                       
            if (K_random > 0)
                index = index((K_nn+2):end, :); % Random Sample from points other than the nearest neighbors but with high affinity
                sorted = sorted((K_nn+2):end, :);
                modval = sum(sorted > thresh);
                modval(modval == 0) = 1;


                % Since we know the number of random neighbors we can
                % preallocate matrices
                %ix_r = repmat(1:n,K_random,1);            
                %rix = bsxfun(@mod, randi(n-(K_nn+1), K_random, n), modval)+1; % Make sure random indices have val > thresh
                %jx_r = index(sub2ind(size(index), rix, ix_r));

                c = 0;
                nn = min(modval, K_random);
                ix_r = zeros(sum(nn),1);
                jx_r = zeros(sum(nn),1);
                for r = 1:n
                    ix_r(c + (1:nn(r))) = r;
                    if (modval(r) == 1)
                        jx_r(c+(1:nn(r))) = index(1,r);
                    else
                        jx_r(c + (1:nn(r))) = randsample(index(1:modval(r),r), nn(r));
                    end
                    c = c+nn(r);
                end

                if (binary)
                    v_r = ones(length(ix_r),1);
                else
                    v_r = W(sub2ind(size(W), jx_r, ix_r));
                end     
            else
                ix_r = zeros(0,1);
                jx_r = zeros(0,1);
                v_r = zeros(0,1);
            end
            ix = [ix_nn(:); ix_r(:)];
            jx = [jx_nn(:); jx_r(:)];
            v = [v_nn(:); v_r(:)];
            Wnn = sparse(jx(:), ix(:), v(:),n,n);
            if (binary)
                Wnn = max(Wnn,Wnn');
            else
                Wnn =.5*(Wnn+Wnn');
            end

        end
        
        function Wnn=compute_W_nn2(W,K, binary)
            [sorted,index] = sort(W,'descend');
            n=length(W);
            
            % Since we know the number of nearest neighbors we can
            % preallocate matrices
            ix = repmat(1:n,K,1);
            jx = index(2:(1+K),:);
            if (binary)
                v = ones(K,n);
            else
                v = sorted(2:(1+K),:);
            end
            
            Wnn = sparse(ix(:), jx(:), v(:),n,n);
            if (binary)
                Wnn = max(Wnn,Wnn');
            else
                Wnn =.5*(Wnn+Wnn');
            end
        end
        
        function WNN=compute_W_nn(W,K)
            [sorted,index] = sort(W,'descend');
            neighborhood = index(2:(1+K),:);
            vals=sorted(2:(1+K),:);
            vals=vals';
            neighborhood=neighborhood';
            D=length(W);
            WNN=sparse(D,D);
            
            for ii=1:D
                jj=neighborhood(ii,:);
                WNN(ii,jj)=vals(ii,:);
            end
            WNN=.5*(WNN+WNN');
        end
        
        function n2 = dist2(x, c)
            %DIST2	Calculates squared distance between two sets of points.
            %
            %	Description
            %	D = DIST2(X, C) takes two matrices of vectors and calculates the
            %	squared Euclidean distance between them.  Both matrices must be of
            %	the same column dimension.  If X has M rows and N columns, and C has
            %	L rows and N columns, then the result has M rows and L columns.  The
            %	I, Jth entry is the  squared distance from the Ith row of X to the
            %	Jth row of C.
            %
            %	See also
            %	GMMACTIV, KMEANS, RBFFWD
            %
            
            %	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)
            
            [ndata, dimx] = size(x);
            [ncentres, dimc] = size(c);
            if dimx ~= dimc
                error('Data dimension does not match dimension of centres')
            end
            
            n2 = (ones(ncentres, 1) * sum((x.^2)', 1))' + ...
                ones(ndata, 1) * sum((c.^2)',1) - ...
                2.*(x*(c'));
            
            % this modification is added to avoid generating negative numbers because of numerical errors.
            
            n2(n2<0)=0;
        end
        
        function fv = flowvar(flow, winsize)
            fv = zeros(size(flow));
            
            % Compute different convolutions
            H = ones(winsize);
            
            % Compute \Sum_{j\inW} ones
            count = conv2(ones(size(flow,1), size(flow,2)), H, 'same');
            
            for i=1:2
                % Compute \Sum_{j\inW} flow_j == \Sum_{j\inW} bar{flow}_i
                sumflow = conv2(flow(:,:,i), H, 'same');
                
                
                % Compute bar{flow} : mean of flow in window
                flowbar = sumflow ./ count;
                
                % Compute \Sum_{j\inW} flow_j.^2
                flow2 = conv2(flow(:,:,i).^2, H, 'same');
                
                % Compute flowvar
                fv(:,:,i) = (flow2 - (2* flowbar .* sumflow) + ((sumflow.^2)./count)) ./ (count - 1);
            end
        end
        
        function fv = flowvar2(flow)
            if (~exist('anisodiff2D','builtin'))
                addpath('/home/elqursh/Projects/Research/Libraries/Filters/AnisoDiff/');
            end
            % Parameters
            num_iter = 30;
            delta_t = 1/7;
            kappa = 30;
            option = 2;
            
            fv = zeros(size(flow));
            for i=1:2
                ad = anisodiff2D(flow(:,:,i),num_iter,delta_t,kappa,option).^2;
                ad2 = anisodiff2D(flow(:,:,i).^2,num_iter,delta_t,kappa,option);
                fv(:,:,i) = abs(ad2 - ad);
            end
            assert(min(fv(:)) >= 0);
        end
        
        function x = solveiter(LU, D ,b)
            assert(size(LU,1) == size(LU,2));
            n_iter = 50;
            [m, n] = size(b);
            x_old = zeros(m,n);
            x = zeros(m,n);
            for i =1:n_iter
                for j=1:m
                    x(j,:) = (b(j,:) - LU(j,:) * x_old) ./ D(j,j);
                end
                x_old = x;
            end
        end
        
        function x = solveiter2(LU ,b)
            assert(size(LU,1) == size(LU,2));
            n_iter = 50;
            [m, n] = size(b);
            x = zeros(m,n);
            for i =1:n_iter
                x = (b - LU * x);
            end
        end
        
        
        function [bestt, bestval, best] = bestsplit(d,D, W, v, l)
            % Find the best way to partition the graph using l equally spaced splitting points
            % d,D : Matricies passed in for efficiency
            % W : Weight matrix (sparse)
            % v : Contineous valued partition
            % l : Number of partitions to
            
            a = max(v);  %Supposed to be one
            b = min(v);  %Supposed to be -ve one
            d2 = (a-b)/(l+1);
            
            p = zeros(length(l),1);
            for i=1:l
                t = b + i*d2;
                p(i) = moseg.MosegUtils.ncut(t,v,W,D,d);
            end
            [bestval,best] = min(p);
            bestt = b + best * d2;
        end
        
        function ncut = ncut(t, v, W, D,d)
            x = (v >= t);
            x = (2 * x) - 1;
            %d = diag(D);
            k = sum(d(x > 0)) / sum(d);
            b = k / (1 - k);
            y = (1 + x) - b * (1 - x);
            ncut = (y' * (D - W) * y) / ( y' * D * y );
        end
        
        function [V,d,D,D3] = segment2(W,k)
            mn = size(W,1);
            
            % Calculate d,b
            d = sum(W,2);
            
            % Calculate A,
            % It is easy to find the exp(D,-1/2) by first using sqrt then find the
            % reciprocal of the diagonal elements
            D = spdiags(d,0,mn,mn);
            %D2 = spdiags(ones(mn,1)./ sqrt(d),0,mn,mn);
            %A = D2 * (D-W) * D2;
            
            % Find k's smallest eigen values
            % Use the generalized egien value problem format :)
            [V,D3] = eigs(D-W,D,k,'sm');
            % Flip the columns of V
            V = fliplr(V);
            D3 = fliplr(flipud(D3));
            %[V,D3] = eigs(A,k,'sm');
        end
        
        function flow = readMiddlebury(path)
            fid = fopen(path, 'r');
            fread(fid, 1, 'float');
            sz = fread(fid, 2, 'int'); %sz = [x y] = [c r]
            flow = fread(fid, inf, 'float');
            flow = permute(reshape(flow, [2 sz(1) sz(2)]),[3 2 1]);
            fclose(fid);
        end
    end
end