classdef MosegAlgLblPropOptions < moseg.MosegAlgAffBaseOptions
    properties 
        %%%% 5. Options for label propagation %%%%%%%%
        lblpropglc = false;
        trackprob = false;
        ZeroOutRangeSP = true; % out of this range set weight to 0. 
        Sp_range = 0.05; % out of this range set weight to 0. 
        
        % a) GLC Options
        n_iter=20;%label propag
        alpha2=.05;% label propagation param
        
        % b) Harmonic Options
        fixlbl = false; % Fix initial probabilities or allow them to change.
        initprob = true; % Initialize unlabled probabilities with previous probability.
        eta = 0;
        eta2 = 0.1;
        normalize = false;        
       
    end
    
    methods
        function obj = MosegAlgLblPropOptions
        end
    end    
end