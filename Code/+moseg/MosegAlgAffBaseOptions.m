classdef MosegAlgAffBaseOptions
    properties
        verbose = true;
        lookahead = 5;
        buffSize = 5; % Number of points to hold in memory for each trajectory.
        smoothSize = 3; % Smallest length of a trajectory for it to be added in embedding.        
        initframes = 3; % How many frames before initializing.

        %%%% 3. Options for distance computation. %%%%
        lrgNum = 100000000000000;
        mode = 'spatial+traj'; % Possible values 'spatial', 'spatial+traj', 'traj'
        conflen = 3;

        %%% Sparsification 
        sparsify = false;
        sparsify_Wthresh = 100 * eps; % If affinity is less than this value then sparsify        
        sparsification_method = 'delaunay'; % 'knn' or 'radius', 'delaunay'
        sparsify_radius = 5; % 20 pixel radius (applicable to sparsification_method = 'radius')
        % Relation: thisval = exp(-Dthresh./sigma_traj), Dthresh = -log(thisval)

        %%%% 4. Options for Affinity computation %%%%
        afffunc = @(x)x;
        trajlen_reweight = false;
%         ZeroOutRangeSP = true; % out of this range set weight to 0. 
%         Sp_range = 3; % out of this range set weight to 0. 
        
        % afffunc = @(x) exp(-x);
        sigma_traj = 0.5;
        sigma_noise = 1; % Noise have a variance of simga_noise pixel/frame, will cap minimum sigma value based on smoothSize * thisval 
        sigma_sp = 300; % Was 100 , to fix problem with spatial propagation
        k_nn = 40; % 4 
        k_random = 40;     % NEW: Select 40 more other the nn in random
        nn_thresh = 0.01;  % NEW
        binW = false;
 
        ncut = true;
        ncut_thresh = 0.0001;        
        outlier_size = 5;
    end
    
    methods
        function obj = MosegAlgAffBaseOptions
        end
    end    
end