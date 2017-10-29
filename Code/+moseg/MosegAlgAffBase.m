% Factor common code between CRF and LblProp here
classdef MosegAlgAffBase < moseg.MosegAlg
    properties
        PrevState;
        CurrState;        
        InitAlg = moseg.MosegAlgDummy;
        Initialized = false;
        StartFrame;
        Word;                
    end

    events
        SelectTrajPair    
    end
        
    methods
        function obj = MosegAlgAffBase(startframe, varargin)            
            obj.SupportsPartialTracks = true;
            
            if (mod(length(varargin),2) ~= 0)
                error('Invalid Number of arguments');
            end
            
            for i=1:2:length(varargin)
                switch varargin{i}
                    case 'InitAlg'
                        obj.InitAlg = varargin{i+1};
                    case 'Options'
                        obj.Options = varargin{i+1};
                end
            end
            
            obj.reset(startframe);            
        end

        function reset(obj, startframe)
            obj.Initialized = false;
            obj.StartFrame = startframe;
            
            obj.CurrState = moseg.MosegAlgAffBaseState(startframe - 1);
            obj.CurrState.frameNo = startframe - 1; 
            obj.PrevState = [];
        end

    end
    
    methods(Access= protected, Abstract=true)
        stepProp(obj)
    end
%     
    methods(Access = protected)                
        function stepInit(obj)
            obj.PrevState = obj.CurrState;            
            obj.CurrState = moseg.MosegAlgAffBaseState(obj.PrevState.frameNo + 1);
            obj.PrevState.W = []; % Clear Memory
        end        
        
        function lbls = stepPartial(obj, points, trajIds, trackno, len, ~)
            obj.TimingStats = struct;
            obj.stepInit();
            tic;
            obj.stepTraj(trajIds, trackno, len);
            obj.TimingStats.traj = toc;

            tic;
            if (obj.Options.sparsify)
                obj.stepDsparse(points, trajIds);
            else
                obj.stepD(points, trajIds);
            end
            obj.TimingStats.dist = toc;
            
            tid2=tic;
            % Compute Affinity matrix
            obj.stepAff();
            obj.TimingStats.aff = toc(tid2);
                 
            lbls = zeros(1, size(points,2));                        

            
            if (~obj.Initialized && (obj.CurrState.frameNo - obj.StartFrame + 1) >= obj.Options.initframes)
                lbls1 = obj.InitAlg.step(points, trajIds, trackno, len, obj.CurrState.frameNo - obj.StartFrame + 1);                
                trajIds = trajIds(lbls1 > 0);
                lbls2 = lbls1(lbls1 > 0);
                
                % Assign labels for points
                [~, ia, ib] = intersect(trajIds, obj.CurrState.memTrajIds);
                obj.CurrState.memLbls(ib) = lbls2(ia) + 1; % Reserve the label 1 for outliers
                
                % Assign probabilities
                m = max(lbls2)+1;
                nl = length(ib);
                obj.CurrState.memProb = zeros(length(obj.CurrState.memTrajIds),m);
                obj.CurrState.memProb(ib,:) = zeros(nl,m);
                obj.CurrState.memProb(sub2ind([nl m], ib', lbls2(ia)+1)) = 1;

                obj.Initialized = true;
            end
            
            if (obj.Initialized) 
                tid = tic;
                obj.stepProp;
                if (obj.Options.ncut == true)            
                    obj.stepNcut;
                end
                obj.TimingStats.prop = toc(tid);

                % Assign labels for points
                [~, ia, ib] = intersect(trajIds, obj.CurrState.memTrajIds);
                lbls(ia) = obj.CurrState.memLbls(ib);

            end

        end

        function stepTraj(obj, trajIds, trackno, len)
            % Version that holds start and ending frames for memory trajectories.
            
            newMemTraj = trajIds(len == obj.Options.smoothSize);
            newMemTrackno = trackno(len == obj.Options.smoothSize);
            %newLbls = obj.CurrState.lbls(obj.CurrState.len == obj.Options.smoothSize);
            newLbls = zeros(1, length(newMemTraj));
            obj.CurrState.memTrackno = [obj.PrevState.memTrackno newMemTrackno];
            obj.CurrState.memTrajIds = [obj.PrevState.memTrajIds newMemTraj];
            obj.CurrState.memLbls = [obj.PrevState.memLbls newLbls];
            m = size(obj.PrevState.memProb,2);
            newProb = zeros(length(newMemTraj),m);                        
            obj.CurrState.memProb = [obj.PrevState.memProb; newProb];
            obj.CurrState.memStartFrame = [obj.PrevState.memStartFrame repmat(obj.CurrState.frameNo - obj.Options.smoothSize +1, 1, length(newMemTraj))];
            obj.CurrState.memEndFrame   = [obj.PrevState.memEndFrame   repmat(obj.CurrState.frameNo, 1, length(newMemTraj))];

            [~,ia,~] = intersect(obj.CurrState.memTrajIds, trajIds);
            obj.CurrState.memEndFrame(ia) = obj.CurrState.frameNo;

            % Drop trajectories ending before current frame,
            keepIdx = (obj.CurrState.memEndFrame == obj.CurrState.frameNo);
            obj.CurrState.memLbls = obj.CurrState.memLbls(keepIdx);
            obj.CurrState.memProb = obj.CurrState.memProb(keepIdx,:);
            obj.CurrState.memTrajIds = obj.CurrState.memTrajIds(keepIdx);
            obj.CurrState.memStartFrame = obj.CurrState.memStartFrame(keepIdx);
            obj.CurrState.memEndFrame = obj.CurrState.memEndFrame(keepIdx);
            obj.CurrState.memTrackno = obj.CurrState.memTrackno(keepIdx);            
        end

        function stepDsparse(obj, points, trajIds)
            % Algorithm
            % Step 1: Precompute motion vectors and points
            % Step 2: Initialize data structures using previous state data
            % Step 3: Update Dm, Dsp matrices only for 1 entries in the sparse
            % mask.
            % Step 4: Compute affinities
            
            n = length(obj.CurrState.memTrajIds);            
            if (n == 0)
                return;
            end             
            
            % To Simplify: Reorder input points to match the order we
            % maintain locally.
            [~, ia, ib] = intersect(obj.CurrState.memTrajIds, trajIds);
            assert(all(ia == (1:n)));
            points = points(:,ib);
            
            % Compute variance of delta and use it to scale distances
            delta = points(1:2,:) - points((2*obj.Options.smoothSize-1):(2*obj.Options.smoothSize),:);
            stddev = std(delta, 0, 2);
            stddev = max(stddev, obj.Options.smoothSize * repmat( obj.Options.sigma_noise, 2,1));

            % Scale distances so that we use the standard gaussian
            delta(1,:) = delta(1,:) ./ (stddev(1) .* obj.Options.sigma_traj);
            delta(2,:) = delta(2,:) ./ (stddev(2) .* obj.Options.sigma_traj);
            points = points ./ obj.Options.sigma_sp;
            
            [~, ia, ib] = intersect(obj.PrevState.memTrajIds, obj.CurrState.memTrajIds);
            [ix, jx] = find(obj.PrevState.Mask(ia,ia));
            ix = ib(ix); jx = ib(jx); 
            % Create a constrained delaunay traingulation from points in
            % current frame
            obj.CurrState.Mask = logical(sparse(ix', jx', 1, n,n));
            switch(obj.Options.sparsification_method)
                case 'delaunay'
                    dt = DelaunayTri(points(1:2,:)');
                    e = dt.edges();
                case 'knn'
                case 'radius'
                    r = obj.Options.sparsify_radius;
                    idx = rangesearch(points(1:2,:)',points(1:2,:)',r);
                    num = num2cell((1:length(idx))');
                    e1 = cellfun(@(x,y) [repmat(x,length(y),1) y'], num, idx, 'UniformOutput',false);
                    e = cell2mat(e1);
            end
            obj.CurrState.Mask = obj.CurrState.Mask | logical(sparse(e(:,1), e(:,2), 1, n,n));
                        
            [ix, jx, ~] = find(obj.CurrState.Mask);            
            if (strcmp(obj.Options.mode, 'traj') || strcmp(obj.Options.mode, 'spatial+traj'))            
                obj.CurrState.Dm = sparse(n,n);            
                obj.CurrState.Dm(ib,ib) = obj.PrevState.Dm(ia,ia);   
                DeltaDm  = sparse(ix, jx,sum(delta(:,ix).^2  + delta(:,jx).^2  - 2 .* delta(:,ix)  .* delta(:,jx)),n,n);
                obj.CurrState.Dm  = max(obj.CurrState.Dm , DeltaDm);
                clear DeletaDm;
            end
            
            if (strcmp(obj.Options.mode, 'spatial') || strcmp(obj.Options.mode, 'spatial+traj'))            
                obj.CurrState.Dsp = sparse(n,n);            
                obj.CurrState.Dsp(ib,ib) = obj.PrevState.Dsp(ia,ia);
                DeltaDsp = sparse(ix, jx, sum(points(1:2,ix).^2 + points(1:2,jx).^2 - 2 .* points(1:2,ix) .* points(1:2,jx)),n,n);
                obj.CurrState.Dsp = max(obj.CurrState.Dsp, DeltaDsp);
                clear DeletaDsp;
            end
        end
                
        function stepD(obj, points, trajIds)
            % Algorithm
            % Step 1: Precompute motion vectors and points
            % Step 2: Initialize data structures using previous state data
            % Step 3: Update Dm, Dsp matrices only for 1 entries in the sparse
            % mask.
            % Step 4: Compute affinities    
            n = length(obj.CurrState.memTrajIds);            
            if (n == 0)
                return;
            end            
            
            % To Simplify: Reorder input points to match the order we
            % maintain locally.

            [~, ia, ib] = intersect(obj.CurrState.memTrajIds, trajIds);
            assert(all(ia' == (1:n)));
            points = points(:,ib);
            
            % Compute variance of delta and use it to scale distances
            delta = points(1:2,:) - points((2*obj.Options.smoothSize-1):(2*obj.Options.smoothSize),:);
            stddev = std(delta, 0, 2);
            stddev = max(stddev, obj.Options.smoothSize * repmat( obj.Options.sigma_noise, 2,1));

            % Scale distances so that we use the standard gaussian
            delta(1,:) = delta(1,:) ./ (stddev(1) .* obj.Options.sigma_traj);
            delta(2,:) = delta(2,:) ./ (stddev(2) .* obj.Options.sigma_traj);
            points = points ./ obj.Options.sigma_sp;
            
            [~, ia, ib] = intersect(obj.PrevState.memTrajIds, obj.CurrState.memTrajIds);

            if (strcmp(obj.Options.mode, 'traj') || strcmp(obj.Options.mode, 'spatial+traj'))
                obj.CurrState.Dm = zeros(n,n);
                obj.CurrState.Dm(ib,ib) = obj.PrevState.Dm(ia,ia);
                obj.PrevState.Dm = [];  % Clear Memory by emptying previous arrays                
                DeltaDm = moseg.MosegUtils.dist2(delta',delta');
                obj.CurrState.Dm = max(obj.CurrState.Dm, DeltaDm);
                DeltaDm = []; % Clear memory
            end
            
            if (strcmp(obj.Options.mode, 'spatial') || strcmp(obj.Options.mode, 'spatial+traj'))
                obj.CurrState.Dsp = zeros(n,n);            
                obj.CurrState.Dsp(ib,ib) = obj.PrevState.Dsp(ia,ia);
                obj.PrevState.Dsp = []; % Clear Memory by emptying previous arrays      
                DeltaDsp = moseg.MosegUtils.dist2(points(1:2,:)', points(1:2,:)');
                obj.CurrState.Dsp = max((obj.CurrState.Dsp ~= obj.Options.lrgNum) .* obj.CurrState.Dsp,DeltaDsp);
                obj.CurrState.DspInRange = obj.CurrState.Dsp < obj.Options.Sp_range; 
               
            end
        end   
        
        function stepAff(obj)
            n = length(obj.CurrState.memTrajIds);  
            fprintf('memTrajIds: %d\n', n); 
            
            if (n == 0)
                return;
            end                        
            
            if (obj.Options.sparsify)            
                % Compute affinities
                obj.CurrState.W = sparse(n,n);
                switch(obj.Options.mode)
                    case 'spatial'
                        obj.CurrState.W(obj.CurrState.Mask) = ...
                            exp(-obj.CurrState.Dsp(obj.CurrState.Mask));
                    case 'spatial+traj'
                        obj.CurrState.W(obj.CurrState.Mask) = ...
                            exp(-( obj.CurrState.Dsp(obj.CurrState.Mask) + obj.CurrState.Dm(obj.CurrState.Mask)));
                    case 'traj'
                        obj.CurrState.W(obj.CurrState.Mask) = ...
                            exp(-obj.CurrState.Dm(obj.CurrState.Mask));            
                end  

                % Compute new mask and apply it to matrices
                obj.CurrState.Mask = obj.CurrState.W > obj.Options.sparsify_Wthresh;
                obj.CurrState.W   = obj.CurrState.W   .* obj.CurrState.Mask;
                switch(obj.Options.mode)
                    case 'spatial'
                        obj.CurrState.Dsp = obj.CurrState.Dsp .* obj.CurrState.Mask;
                    case 'spatial+traj'
                        obj.CurrState.Dsp = obj.CurrState.Dsp .* obj.CurrState.Mask;                    
                        obj.CurrState.Dm  = obj.CurrState.Dm  .* obj.CurrState.Mask;
                    case 'traj'
                        obj.CurrState.Dm  = obj.CurrState.Dm  .* obj.CurrState.Mask;
                end                  
            else
                switch(obj.Options.mode)
                    case 'spatial'
                        obj.CurrState.W = exp(-obj.CurrState.Dsp);
                    case 'spatial+traj'
                        obj.CurrState.W = exp(-( obj.CurrState.Dsp + obj.CurrState.Dm));                
                    case 'traj'
                        obj.CurrState.W = exp(-obj.CurrState.Dm);                
                end                                
            end
            
            % Reweight affinities based on common length
            if (obj.Options.trajlen_reweight == true)
                n = length(obj.CurrState.memStartFrame);                
                % Weight each entry by the length of overlap
                commonStart = max(repmat(obj.CurrState.memStartFrame',1,n),repmat(obj.CurrState.memStartFrame,n,1));
                commonEnd   = min(repmat(obj.CurrState.memEndFrame',1,n), repmat(obj.CurrState.memEndFrame,n,1));        
                commonLength = max(commonEnd - commonStart,0);    
                % Reweight affinities based on
                obj.CurrState.W = obj.CurrState.W .* commonLength;
            end
            if(obj.Options.ZeroOutRangeSP == true)
                obj.CurrState.W  =  obj.CurrState.W .* obj.CurrState.DspInRange; 
            end
            
            % Sparsify by computing nearest neighbors and zeroing out <
            % 10*eps values
            %obj.CurrState.W = moseg.MosegUtils.compute_W_nn2(obj.CurrState.W, obj.Options.k_nn, obj.Options.binW);
            obj.CurrState.W = moseg.MosegUtils.nnSparseMat(obj.CurrState.W, obj.Options.k_nn, obj.Options.k_random, obj.Options.nn_thresh, obj.Options.binW);
            
            % Handle outliers
            % Zero out small componenets.
            ix = obj.CurrState.W > 0 & obj.CurrState.W < 10 * eps;
            obj.CurrState.W(ix) = 0;            
        end
        
        function stepNcut(obj)
            n = length(obj.CurrState.memTrajIds);

            % See if we need to split any of the clusters
            clusters = unique(obj.CurrState.memLbls);
            clusters = clusters(clusters ~= 0 & clusters ~=1);

            memProb = obj.CurrState.memProb;
            memLbls = obj.CurrState.memLbls;
            numclust = size(memProb,2);
            for i=clusters
                if (sum(memLbls == i) < 15)
%                         fprintf('Segment is too small to split\n');
                    continue;
                end
                ix = obj.evalsplit2(i);
                % Move outliers to the outlier cluster.
                if (~isempty(ix))
                    memProb(ix{1},1) = 1;
                    memProb(ix{1},i) = 0;
                    memLbls(ix{1}) = 1;
                end

                % Cluster 2 retains the same cluster number and the
                % resut are assigned new clusters
%                   fprintf('Splitting cluster %d\n', i);
                for j=3:length(ix)
                    newclustnum = numclust + 1;
                    numclust = numclust + 1;
                    memProb = [memProb zeros(n,1)];
                    memProb(ix{j},newclustnum) = memProb(ix{j}, i);
                    memProb(ix{j}, i) = 0;
                    memLbls(ix{j}) = newclustnum;
                end
            end
            obj.CurrState.memProb = memProb;
            obj.CurrState.memLbls = memLbls;
        end
        
        function ix = evalsplit2(obj, clustnum)
            numccsplit = 0;
            numnormcutsplit = 0;
            numsmall = 0;
            numoutliers = 0;
            
            clustix = find(obj.CurrState.memLbls == clustnum);
            Wc = obj.CurrState.W(clustix, clustix);
            Wceps = Wc > (eps * 10); % Remove v. weak edges to avoid numerical problems
            
            % First find connected components 
            [S,C] = graphconncomp(Wceps);
            % Find the counts in each connected component
            %fprintf('Cluster %d has %d connected componenets\n', clustnum, S);
            ix2 = cell(S,1);
            for i=1:S
                ix2{i} = clustix(C==i);
            end
            numccsplit = numccsplit + S;
            
            ix = cell(2*S,1);
            numclust = 0;
            for i=1:S
                % Ignore small clusters
                if (length(ix2{i}) < 15)
                    ix{numclust+1} = ix2{i};
                    numclust = numclust + 1;
                    numsmall=  numsmall + 1;
                    continue;
                end
                
                clustix2 = ix2{i};
                
                Wc2 = obj.CurrState.W(clustix2, clustix2);
                n = length(clustix2);
                % Solve the normalized eigen value problem
                d = sum(Wc2,2);
                Dm = spdiags(d, 0, n,n);
                L = Dm - Wc2;
                % Since we know that in L has rank n-1 start from a small value
                % of sigma instead of using 0;.
                if ( n < 500)
                    [V, Sigma] = eig(full(L), full(Dm));
                    sig = diag(Sigma);
                else
                    try
                        opts.issym = 1;
                        [V,Sigma] = eigs(L,Dm,2,'sm',opts); % Extract 2 smallest eigen values/vectors
                        sig = spdiags(Sigma,0);            
                    catch err
                        fprintf('Error occured in sparse eigen value decomp.\n');
                        [V, Sigma] = eig(full(L), full(Dm));
                        sig = diag(Sigma);                    
                    end
                end
                [~, isorted] = sort(sig,'ascend');
                v = V(:, isorted(2));            
                [bestt, bestval] = moseg.MosegUtils.bestsplit(d, Dm, Wc2, v, 30);
                split = (bestval < obj.Options.ncut_thresh);
                cix1 = clustix2(v < bestt);
                cix2 = clustix2(v >= bestt);

                if (min(length(cix1), length(cix2)) < 15)
%                     fprintf('Uneven split < 15, not splitting\n');
                    split = false;
                end
                if (split)
                    ix{numclust+1} = cix1;
                    ix{numclust+2} = cix2;
                    numclust = numclust + 2;
                    numnormcutsplit = numnormcutsplit + 1;
                else
                    ix{numclust+1} = clustix2;
                    numclust = numclust +1;
                end
            end

            ix = ix(1:numclust);

            % Find counts 
            counts = cellfun(@length, ix);

            % Give singletons the special label 1 (for outliers)
            numoutliers = sum(counts <= obj.Options.outlier_size);
            outliers = [ix{counts <= obj.Options.outlier_size}];
            nonoutliers = ix(counts > obj.Options.outlier_size);
            
            ix = cell(length(nonoutliers)+1, 1);
            ix(2:end) = nonoutliers;
            ix{1} = outliers;
            %ix = ix(1:sum(counts~=1)+1);

            % Find counts 
            counts = cellfun(@length, ix);
            
            % Swap
            [~,m] = max(counts);
            if (m > 2)
                temp = ix{2};                
                ix{2} = ix{m};
                ix{m} = temp;
            end
            obj.debugmsg('\tClust %d splits: CC = %d, NC = %d, Outlier =%d. Small = %d\n', clustnum, numccsplit, numnormcutsplit, numoutliers, numsmall);
        end        
        
    end
    
end


