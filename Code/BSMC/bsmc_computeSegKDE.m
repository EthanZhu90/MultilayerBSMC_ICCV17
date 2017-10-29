function [state] = bsmc_computeSegKDE(state, options)
    debug_info = struct;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize with Maximum Liklihood segmentation using graph cut.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    layer_num = length(state.layers); 
    label_mask = bsmc_computeGCInitSeg_multilabels(state, options, options.lambda_gc);
    
    %%% this counts is only for segmentation. final count is update later,
    %%% not this one. 
    app_addCurr_count = cell(1, layer_num);
    for i = 1: layer_num
        app_addCurr_count{i} = state.layers(i).app_model.counts + double(label_mask == state.layers(i).label); 
    end

    [nrow, ncol] = size(label_mask);  
    debug_info.InitSegIm = label_mask;
    app_counts = cell(1,layer_num); 
    app_invalid_LsPixel = cell(1, layer_num);  %% for part that 0 < app count < options.kde_start_eval 0 -6
    app_invalid_NoPixel = cell(1, layer_num); 
    for i = 1:layer_num
        app_counts{i} = min(state.layers(i).app_model.counts, options.kde_n); %10
        app_invalid_NoPixel{i} = app_addCurr_count{i} == 0 ; %5
        app_invalid_LsPixel{i} = (app_addCurr_count{i} > 0) & (app_addCurr_count{i} < options.kde_start_eval+1); 
    end
    debug_info.app_invalid_NoPixel = app_invalid_NoPixel;
    debug_info.app_invalid_LsPixel = app_invalid_LsPixel;
    

    % given the predicted label prior and appearance models, a MAP
    % estimate of the labels is inferred.
    for i = 1:layer_num
        state.layers(i).prob_model.post = bsmc_evalKDE_padded(state.layers(i).app_model.x, state.layers(i).app_model.counts,...
            state.obs, state.layers(i).motion_model.var, options);
        state.layers(i).prob_model.post(app_invalid_NoPixel{i}) = options.kde_thresh;  %%% reset prob: step 1
    end

    %%%%%%%% add equal prob mask %%%%% reset prob: step 2
    zeroMask = zeros(nrow, ncol); 
    for i = 1:layer_num
        zeroMask = zeroMask | app_invalid_LsPixel{i};
    end
    for i = 1:layer_num
        state.layers(i).prob_model.post(zeroMask) = 1; 
    end
    debug_info.app_invalid_LsPixel_Com = zeroMask; 

    %%%%% reset prob: step 3: protect the new layers 
    idx = find([state.layers.nframe] <  (options.kde_start_eval+1)); 
    if(~isempty(idx))
        newlayer_mask = zeros(nrow,ncol); 
        for i = idx
            newlayer_mask = newlayer_mask | (label_mask == state.layers(i).label); 
        end
        debug_info.newlayer_mask = newlayer_mask; 
        oneMask = ones(nrow, ncol); 
        for i = 1:layer_num
            state.layers(i).prob_model.post((label_mask == state.layers(i).label) & newlayer_mask) =...
                (options.kde_thresh + options.kde_shift).*oneMask((label_mask == state.layers(i).label) & newlayer_mask); 
            state.layers(i).prob_model.post((label_mask ~= state.layers(i).label) & newlayer_mask) =...
                (options.kde_thresh - options.kde_shift).*oneMask((label_mask ~= state.layers(i).label) & newlayer_mask); 
        end
    end

   
    for i = 1:layer_num
        if i == 1
            post_norm_scale = state.layers(1).prob_model.post; 
        else
            post_norm_scale = post_norm_scale + state.layers(i).prob_model.post;
        end
    end 
    if (options.label_prior == 0)     
        for i = 1:layer_num
            state.layers(i).prob_model.post = state.layers(i).prob_model.post ./ post_norm_scale; 
        end 
    else
        bgprior = (1-state.model.fgprior);
        fgprior = state.model.fgprior;
        bgpost = (bgprior .* bgprob) ./ (bgprior .* bgprob + fgprior .* fgprob);
        
        % For places where we have only one appearance model, likelihood
        % can be very large which gives us one or zeros
        % This prevents the code from recovering at a later iteration.
        % Scale the values to an acceptable range, say between 0.05 and
        % 0.95. Basically here we decide how much temporal smoothing do we
        % want, beside of course the place where we predict the labels.
        bgpost2 = bgprob ./ ( bgprob + fgprob ); 
        bgpost2 = 0.05 + 0.5 .* bgpost2; 
        fgpost2 = 1 - bgpost2;
        bgpost3 = (bgprior .* bgpost2) ./ (bgprior .* bgpost2 + fgprior .* fgpost2);
        a = xor(fginvalid,bginvalid);
        bgpost(a == 1) = bgpost3(a==1);
    end
   
    [m,n,~] = size(state.frame);
    
    post_matrix = zeros(m,n,layer_num);
    for i = 1:layer_num
         post_matrix(:,:,i) = state.layers(i).prob_model.post; 
    end 
    
    [~, mask] = max(post_matrix,[],3); % find the max along dim 3 
    % replace mask with cluster label 
    lmask = zeros(size(mask)); 
    for i = 1:layer_num
        lmask(mask==i) = state.layers(i).label; 
    end 

    state.beforeGC = lmask;  
    if (options.postGC == 1)
        state = computePostGraphCutSegmentation_multilabels(state, lmask, options); %%% label as the data cost
    end
    
    %figure(2)
    %imshow(mask{1}); 
    if (options.postGC == 2)
        state = computePostGraphCutSegmentationProb_multilabels(state, options); %%% prob as the data cost 
    end
    %figure(3)
    %imshow(mask{1});
    state.debug_info = debug_info; 
end


function [state] = computePostGraphCutSegmentationProb_multilabels(state, options) %by ethan
    dframe = double(state.frame);
    [m n ~] = size(dframe); 
    

    layer_num = length(state.layers); 
    %idx = find([state.layers.label] == bgClust); % it's bg.
    initial_labels =randi([1, layer_num],1,m*n);  % for CLASS   
    
    
    %initial_labels = reshape(mask, [1,m*n]);  % for CLASS   
    initial_labels = initial_labels -1; % became 0,1,2
    
    
    %%%%%%%% prepare the data cost (UNARY)
   
 
    layer_num = length(state.layers); 
    datacost = []; 
    for i = 1:layer_num
        post = 1 - state.layers(i).prob_model.post; 
        post = reshape(post, 1, m*n); 
        datacost = [datacost; post]; 
        
    end 
    %%% no need to normal; 
%     norm_scale = sum(datacost,1); 
%     norm_scale = repmat(norm_scale, layer_num, 1); 
%     datacost= datacost./norm_scale; %%% normalization 

    %%%%%%%% prepare label cost (LABELCOST)
    labelcost = ones(layer_num); 
    for i = 1:layer_num
        labelcost(i, i) = 0; 
    end
    
    %%%%%%%% prepare edge weigh cost (PAIRWISE)
    hDiff = [dframe(:,1:n-1,:) - dframe(:,2:n,:), zeros(m,1,3)];
    vDiff = [dframe(1:m-1,:,:) - dframe(2:m,:,:); zeros(1,n,3)];
    seDiff = [dframe(1:m-1,1:n-1,:) - dframe(2:m,2:n,:) zeros(m-1,1,3); zeros(1,n,3)];
    neDiff = [zeros(1,n,3); dframe(2:m,1:n-1,:) - dframe(1:m-1,2:n,:) zeros(m-1,1,3)];
    
    halfInvSig = (1/(options.sigma_app)) * eye(3);
    hZ = halfInvSig * reshape(permute(hDiff, [3 1 2]), 3, m*n);
    vZ = halfInvSig * reshape(permute(vDiff, [3 1 2]), 3, m*n);
    seZ = halfInvSig * reshape(permute(seDiff, [3 1 2]), 3, m*n);
    neZ = halfInvSig * reshape(permute(neDiff, [3 1 2]), 3, m*n);
 
    hCue = reshape(exp(-0.5 * sum(hZ .* hZ)), m, n);
    vCue = reshape(exp(-0.5 * sum(vZ .* vZ)), m, n);
    seCue = reshape(exp(-0.5 * sum(seZ .* seZ)), m,n);  
    neCue = reshape(exp(-0.5 * sum(neZ .* neZ)), m,n);
   
    
    rows_store = []; 
    cols_store = [];
    v_store = []; 
    % prepare pairwise for hCue data 
    hCue = hCue(:, 1:n-1);% m * n-1
    i = 1:m; 
    i = repmat(i, n-1, 1); 
    i = reshape(i, [1, (n-1)*m]); 
    
    j = (1:n-1)';
    j = repmat(j, 1, m); 
    j = reshape(j, [1, (n-1)*m]); 
    
    rows = sub2ind([m,n], i, j);
    cols = sub2ind([m,n], i, j+1); 
    v = reshape(hCue', [1, m*(n-1)]); 
    rows_store = [rows_store, rows];
    cols_store = [cols_store, cols];
    v_store = [v_store, v]; 
    rows_store = [rows_store, cols];
    cols_store = [cols_store, rows];
    v_store = [v_store, v]; 
    
    % prepare pairwise for vCue data 
    vCue = vCue(1:m-1, :); % m-1 * n
    i = 1:m-1; 
    i = repmat(i, n, 1); 
    i = reshape(i, [1, n*(m-1)]); 
    
    j = (1:n)';
    j= repmat(j, 1, m-1);
    j = reshape(j, [1, n*(m-1)]); 
    
    rows = sub2ind([m,n], i, j);
    cols = sub2ind([m,n], i+1, j); 
    v = reshape(vCue', [1, n*(m-1)]); 
    rows_store = [rows_store, rows];
    cols_store = [cols_store, cols];
    v_store = [v_store, v]; 
    rows_store = [rows_store, cols];
    cols_store = [cols_store, rows];
    v_store = [v_store, v]; 
    
    % prepare pairwise for seCue data 
    seCue = seCue(1:m-1, 1:n-1); % m-1 * n-1
    i = 1:m-1; 
    i = repmat(i, n-1, 1); 
    i = reshape(i, [1, (n-1)*(m-1)]); 
    
    j = (1:n-1)';
    j= repmat(j, 1, m-1);
    j = reshape(j, [1, (n-1)*(m-1)]); 
    rows = sub2ind([m,n], i, j);
    cols = sub2ind([m,n], i+1, j+1); 
    v = reshape(seCue', [1, (n-1)*(m-1)]); 
    rows_store = [rows_store, rows];
    cols_store = [cols_store, cols];
    v_store = [v_store, v]; 
    rows_store = [rows_store, cols];
    cols_store = [cols_store, rows];
    v_store = [v_store, v]; 
    
    % prepare pairwise for seCue data 
    neCue = neCue(2:m, 1:n-1); % m-1 * n-1
    i = 2:m; 
    i = repmat(i, n-1, 1); 
    i = reshape(i, [1, (n-1)*(m-1)]); 
    
    j = (1:n-1)';
    j= repmat(j, 1, m-1);
    j = reshape(j, [1, (n-1)*(m-1)]); 
    
    rows = sub2ind([m,n], i, j);
    cols = sub2ind([m,n], i-1, j+1); 
    v = reshape(neCue', [1, (n-1)*(m-1)]); 
    rows_store = [rows_store, rows];
    cols_store = [cols_store, cols];
    v_store = [v_store, v]; 
    rows_store = [rows_store, cols];
    cols_store = [cols_store, rows];
    v_store = [v_store, v];  
    
    sMatrix = sparse(rows_store, cols_store, v_store, m*n, m*n);  % N X N sparse matrix  
      
    %EXPANSION = 0; %  0 == swap(2 lablels), 1 == expansion
    if(layer_num > 2)
        EXPANSION = 1; %  0 == swap(2 lablels), 1 == expansion
    else
        EXPANSION = 0;
    end
    dinitial_labels = double(initial_labels);
    sdatacost = single(datacost);
    slabelcost = single(labelcost); 

    % use graph cut with multiple labels
    
    [labels energy energy_after] = GCMex(dinitial_labels, options.datacost_lambda *sdatacost, ...
        options.pairwise_lambda * sMatrix, options.labelcost_lambda * slabelcost, EXPANSION);
  
    mask = reshape(labels, [m , n]);
    mask = mask + 1; 
   
    label_mask = zeros(size(mask));     
    for i = 1: layer_num
        %mask_term = zeros(size(mask)); 
        %mask_term(mask==i) = 1; 
        %state.layers(i).mask = mask_term; 
        label_mask(mask==i) = state.layers(i).label; 
    end
    state.lmask = label_mask; 
end


function [state] = computePostGraphCutSegmentation_multilabels(state, lmask, options) %by ethan
    dframe = double(state.frame);
    [m n ~] = size(dframe); 
    

    
    layer_num = length(state.layers); 
    %idx = find([state.layers.label] == bgClust); % it's bg.
    initial_labels =randi([1, layer_num],1,m*n);  % for CLASS   
    
    
    %initial_labels = reshape(mask, [1,m*n]);  % for CLASS   
    initial_labels = initial_labels -1; % became 0,1,2
    
    
    %%%%%%%% prepare the data cost (UNARY)
    bmask = zeros(size(lmask)); 
    layer_num = length(state.layers); 
    for i = 1:layer_num
        bmask(lmask==state.layers(i).label) = i; 
    end 
    bmask = reshape(bmask, [1,m*n]); 
    datacost = []; % (C X N) (1,0,1)
    for i = 1: layer_num
        datacost_term = 0.5 * ones(1, m*n);
        datacost_term(bmask == i) = 0.0001; 
        datacost = [datacost; datacost_term]; 
    end
    
    norm_scale = sum(datacost,1); 
    norm_scale = repmat(norm_scale, layer_num, 1); 
    datacost= datacost./norm_scale; %%% normalization 

    %%%%%%%% prepare label cost (LABELCOST)
    labelcost = ones(layer_num); 
    for i = 1:layer_num
        labelcost(i, i) = 0; 
    end
    
    %%%%%%%% prepare edge weigh cost (PAIRWISE)
    hDiff = [dframe(:,1:n-1,:) - dframe(:,2:n,:), zeros(m,1,3)];
    vDiff = [dframe(1:m-1,:,:) - dframe(2:m,:,:); zeros(1,n,3)];
    seDiff = [dframe(1:m-1,1:n-1,:) - dframe(2:m,2:n,:) zeros(m-1,1,3); zeros(1,n,3)];
    neDiff = [zeros(1,n,3); dframe(2:m,1:n-1,:) - dframe(1:m-1,2:n,:) zeros(m-1,1,3)];
    
    halfInvSig = (1/(options.sigma_app)) * eye(3);
    hZ = halfInvSig * reshape(permute(hDiff, [3 1 2]), 3, m*n);
    vZ = halfInvSig * reshape(permute(vDiff, [3 1 2]), 3, m*n);
    seZ = halfInvSig * reshape(permute(seDiff, [3 1 2]), 3, m*n);
    neZ = halfInvSig * reshape(permute(neDiff, [3 1 2]), 3, m*n);
 
    hCue = reshape(exp(-0.5 * sum(hZ .* hZ)), m, n);
    vCue = reshape(exp(-0.5 * sum(vZ .* vZ)), m, n);
    seCue = reshape(exp(-0.5 * sum(seZ .* seZ)), m,n);  
    neCue = reshape(exp(-0.5 * sum(neZ .* neZ)), m,n);
   
    
    rows_store = []; 
    cols_store = [];
    v_store = []; 
    % prepare pairwise for hCue data 
    hCue = hCue(:, 1:n-1);% m * n-1
    i = 1:m; 
    i = repmat(i, n-1, 1); 
    i = reshape(i, [1, (n-1)*m]); 
    
    j = (1:n-1)';
    j = repmat(j, 1, m); 
    j = reshape(j, [1, (n-1)*m]); 
    
    rows = sub2ind([m,n], i, j);
    cols = sub2ind([m,n], i, j+1); 
    v = reshape(hCue', [1, m*(n-1)]); 
    rows_store = [rows_store, rows];
    cols_store = [cols_store, cols];
    v_store = [v_store, v]; 
    rows_store = [rows_store, cols];
    cols_store = [cols_store, rows];
    v_store = [v_store, v]; 
    
    % prepare pairwise for vCue data 
    vCue = vCue(1:m-1, :); % m-1 * n
    i = 1:m-1; 
    i = repmat(i, n, 1); 
    i = reshape(i, [1, n*(m-1)]); 
    
    j = (1:n)';
    j= repmat(j, 1, m-1);
    j = reshape(j, [1, n*(m-1)]); 
    
    rows = sub2ind([m,n], i, j);
    cols = sub2ind([m,n], i+1, j); 
    v = reshape(vCue', [1, n*(m-1)]); 
    rows_store = [rows_store, rows];
    cols_store = [cols_store, cols];
    v_store = [v_store, v]; 
    rows_store = [rows_store, cols];
    cols_store = [cols_store, rows];
    v_store = [v_store, v]; 
    
    % prepare pairwise for seCue data 
    seCue = seCue(1:m-1, 1:n-1); % m-1 * n-1
    i = 1:m-1; 
    i = repmat(i, n-1, 1); 
    i = reshape(i, [1, (n-1)*(m-1)]); 
    
    j = (1:n-1)';
    j= repmat(j, 1, m-1);
    j = reshape(j, [1, (n-1)*(m-1)]); 
    rows = sub2ind([m,n], i, j);
    cols = sub2ind([m,n], i+1, j+1); 
    v = reshape(seCue', [1, (n-1)*(m-1)]); 
    rows_store = [rows_store, rows];
    cols_store = [cols_store, cols];
    v_store = [v_store, v]; 
    rows_store = [rows_store, cols];
    cols_store = [cols_store, rows];
    v_store = [v_store, v]; 
    
    % prepare pairwise for seCue data 
    neCue = neCue(2:m, 1:n-1); % m-1 * n-1
    i = 2:m; 
    i = repmat(i, n-1, 1); 
    i = reshape(i, [1, (n-1)*(m-1)]); 
    
    j = (1:n-1)';
    j= repmat(j, 1, m-1);
    j = reshape(j, [1, (n-1)*(m-1)]); 
    
    rows = sub2ind([m,n], i, j);
    cols = sub2ind([m,n], i-1, j+1); 
    v = reshape(neCue', [1, (n-1)*(m-1)]); 
    rows_store = [rows_store, rows];
    cols_store = [cols_store, cols];
    v_store = [v_store, v]; 
    rows_store = [rows_store, cols];
    cols_store = [cols_store, rows];
    v_store = [v_store, v];  
    
    sMatrix = sparse(rows_store, cols_store, v_store, m*n, m*n);  % N X N sparse matrix  
      
    %EXPANSION = 0; %  0 == swap(2 lablels), 1 == expansion
    if(layer_num > 2)
        EXPANSION = 1; %  0 == swap(2 lablels), 1 == expansion
    else
        EXPANSION = 0;
    end
    dinitial_labels = double(initial_labels);
    sdatacost = single(datacost);
    slabelcost = single(labelcost); 

    % use graph cut with multiple labels
    
    [labels energy energy_after] = GCMex(dinitial_labels, options.datacost_lambda *sdatacost, ...
        options.pairwise_lambda * sMatrix, options.labelcost_lambda * slabelcost, EXPANSION);
  
    mask = reshape(labels, [m , n]);
    mask = mask + 1; 
   
    label_mask = zeros(size(mask));     
    for i = 1: layer_num
        %mask_term = zeros(size(mask)); 
        %mask_term(mask==i) = 1; 
        %state.layers(i).mask = mask_term; 
        label_mask(mask==i) = state.layers(i).label; 
    end
    state.lmask = label_mask; 
end

function mask = computePostGraphCutSegmentationProb(state, bgpost, options)
    dframe = double(state.frame);
    [m n ~] = size(dframe);            
    hDiff = [dframe(:,1:n-1,:) - dframe(:,2:n,:) zeros(m,1,3)];
    vDiff = [dframe(1:m-1,:,:) - dframe(2:m,:,:); zeros(1,n,3)];
    seDiff = [dframe(1:m-1,1:n-1,:) - dframe(2:m,2:n,:) zeros(m-1,1,3); zeros(1,n,3)];
    neDiff = [zeros(1,n,3); dframe(2:m,1:n-1,:) - dframe(1:m-1,2:n,:) zeros(m-1,1,3)];
    halfInvSig = (1/(options.sigma_app)) * eye(3);
    hZ = halfInvSig * reshape(permute(hDiff, [3 1 2]), 3, m*n);
    vZ = halfInvSig * reshape(permute(vDiff, [3 1 2]), 3, m*n);
    neZ = halfInvSig * reshape(permute(neDiff, [3 1 2]), 3, m*n);
    seZ = halfInvSig * reshape(permute(seDiff, [3 1 2]), 3, m*n);
    hCue = reshape(exp(-0.5 * sum(hZ .* hZ)), m, n);
    vCue = reshape(exp(-0.5 * sum(vZ .* vZ)), m, n);
    neCue = reshape(exp(-0.5 * sum(neZ .* neZ)), m,n);
    seCue = reshape(exp(-0.5 * sum(seZ .* seZ)), m,n);            
    dataCost(:,:,1) = 1-bgpost;
    dataCost(:,:,2) = bgpost;
    mask = bsmc_gridCut(2 * dataCost,  8 * hCue, 8 *vCue, 8 *seCue, 8 *neCue);
end


function mask = computePostGraphCutSegmentation(state, mask, options)
    dframe = double(state.frame);
    [m n ~] = size(dframe);            
    hDiff = [dframe(:,1:n-1,:) - dframe(:,2:n,:) zeros(m,1,3)];
    vDiff = [dframe(1:m-1,:,:) - dframe(2:m,:,:); zeros(1,n,3)];
    seDiff = [dframe(1:m-1,1:n-1,:) - dframe(2:m,2:n,:) zeros(m-1,1,3); zeros(1,n,3)];
    neDiff = [zeros(1,n,3); dframe(2:m,1:n-1,:) - dframe(1:m-1,2:n,:) zeros(m-1,1,3)];
    halfInvSig = (1/(options.sigma_app)) * eye(3);
    hZ = halfInvSig * reshape(permute(hDiff, [3 1 2]), 3, m*n);
    vZ = halfInvSig * reshape(permute(vDiff, [3 1 2]), 3, m*n);
    neZ = halfInvSig * reshape(permute(neDiff, [3 1 2]), 3, m*n);
    seZ = halfInvSig * reshape(permute(seDiff, [3 1 2]), 3, m*n);
    hCue = reshape(exp(-0.5 * sum(hZ .* hZ)), m, n);
    vCue = reshape(exp(-0.5 * sum(vZ .* vZ)), m, n);
    neCue = reshape(exp(-0.5 * sum(neZ .* neZ)), m,n);
    seCue = reshape(exp(-0.5 * sum(seZ .* seZ)), m,n);     

    sparse_fgll_noapp = 0.5 * ones(size(state.frame,1), size(state.frame,2));
    sparse_bgll_noapp = 0.5 * ones(size(state.frame,1), size(state.frame,2));
    
    sparse_fgll_noapp(mask) = options.lrgNum;
    sparse_bgll_noapp(~mask) = options.lrgNum;
    norm_factor = sparse_fgll_noapp + sparse_bgll_noapp;
    sparse_fgll_noapp = sparse_fgll_noapp ./ norm_factor;
    sparse_bgll_noapp = sparse_bgll_noapp ./ norm_factor;
    
    dataCost(:,:,1) = sparse_fgll_noapp;
    dataCost(:,:,2) = sparse_bgll_noapp;
    mask = bsmc_gridCut(2 * dataCost,  8 * hCue, 8 *vCue, 8 *seCue, 8 *neCue);
end

function mask = computeGraphCutSegmentation(state,options)
    dframe = double(state.frame);
    [m n ~] = size(dframe);            
    hDiff = [dframe(:,1:n-1,:) - dframe(:,2:n,:) zeros(m,1,3)];
    vDiff = [dframe(1:m-1,:,:) - dframe(2:m,:,:); zeros(1,n,3)];
    seDiff = [dframe(1:m-1,1:n-1,:) - dframe(2:m,2:n,:) zeros(m-1,1,3); zeros(1,n,3)];
    neDiff = [zeros(1,n,3); dframe(2:m,1:n-1,:) - dframe(1:m-1,2:n,:) zeros(m-1,1,3)];
    halfInvSig = (1/(options.sigma_app)) * eye(3);
    hZ = halfInvSig * reshape(permute(hDiff, [3 1 2]), 3, m*n);
    vZ = halfInvSig * reshape(permute(vDiff, [3 1 2]), 3, m*n);
    neZ = halfInvSig * reshape(permute(neDiff, [3 1 2]), 3, m*n);
    seZ = halfInvSig * reshape(permute(seDiff, [3 1 2]), 3, m*n);
    hCue = reshape(exp(-0.5 * sum(hZ .* hZ)), m, n);
    vCue = reshape(exp(-0.5 * sum(vZ .* vZ)), m, n);
    neCue = reshape(exp(-0.5 * sum(neZ .* neZ)), m,n);
    seCue = reshape(exp(-0.5 * sum(seZ .* seZ)), m,n);            

%     [~, ia,ib] = intersect(state.memTrajIds, state.trajIds);
%     sparse_fgpoints = round(state.points(1:2, ib(state.clust.lbls(ia) == state.clust.fgClust)));
%     sparse_bgpoints = round(state.points(1:2, ib(state.clust.lbls(ia) ~= state.clust.fgClust)));
    sparse_fgpoints = state.sparse_fgpoints;
    sparse_bgpoints = state.sparse_bgpoints;

    sparse_fgll_noapp = 0.5 * ones(size(state.frame,1), size(state.frame,2));
    sparse_bgll_noapp = 0.5 * ones(size(state.frame,1), size(state.frame,2));
    sparse_fgll_noapp(sub2ind(size(sparse_fgll_noapp), sparse_fgpoints(2,:), sparse_fgpoints(1,:))) = options.lrgNum;
    sparse_bgll_noapp(sub2ind(size(sparse_bgll_noapp), sparse_bgpoints(2,:), sparse_bgpoints(1,:))) = options.lrgNum;
    norm_factor = sparse_fgll_noapp + sparse_bgll_noapp;
    seg.sparse_fgll_noapp = sparse_fgll_noapp ./ norm_factor;
    seg.sparse_bgll_noapp = sparse_bgll_noapp ./ norm_factor;
    dataCost(:,:,1) = seg.sparse_fgll_noapp;
    dataCost(:,:,2) = seg.sparse_bgll_noapp;
    mask = bsmc_gridCut(4 * dataCost,  8 * hCue, 8 *vCue, 8 *seCue, 8 *neCue);
end