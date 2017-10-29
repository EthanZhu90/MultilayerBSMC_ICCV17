function label_mask = bsmc_computeGCInitSeg_multilabels( state, options, lambda) %by ethan
%BSMC_COMPUTEGCINITSEG_MULTILABELS Summary of this function goes here
%   Detailed explanation goes here
    addpath('../Utils/GraphCutMex/gco-v2.3');
    bgClust = state.bgClust;
    dframe = double(state.frame);
    [m, n, ~] = size(dframe);
    
    layer_num = length(state.layers); 
    %idx = find([state.layers.label] == bgClust); % it's bg.
    initial_labels =randi([1, layer_num],1,m*n);  % for CLASS   
    
    
    layers = state.layers; 
    for i = 1: layer_num
      %  if(layers(i).label ~= bgClust) % if it's fg.
            sparse_points = layers(i).traj_store.sparse_points; 
            initial_labels(sub2ind([m,n], sparse_points(2,:), sparse_points(1,:))) = i;   
      %  end
    end
    initial_labels = initial_labels -1; 
    %%%%%%%% prepare the data cost (UNARY)
    datacost = 0.5 * ones(layer_num, m*n); % (C X N) (1,0,1)
    for i = 1: layer_num
        sparse_points = layers(i).traj_store.sparse_points;
        traj_num = size(sparse_points,2); 
        datacost_term = 0.5 * ones(layer_num, 1); 
        datacost_term(i) = 0.0001; %this is a right way %options.lrgNum; 
        datacost_term = repmat(datacost_term, 1, traj_num); 
        datacost(:, sub2ind([m,n], sparse_points(2,:), sparse_points(1,:))) = datacost_term;  
    end
    norm_scale = sum(datacost,1); 
    norm_scale = repmat(norm_scale, layer_num, 1); 
    datacost= datacost./norm_scale; 
    
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
 
   [labels energy energy_after] = GCMex(dinitial_labels, options.InGC_datacost_lambda *sdatacost, ...
        options.InGC_pairwise_lambda * sMatrix, options.InGC_labelcost_lambda * slabelcost, EXPANSION);
   mask = reshape(labels, [m , n]);
   mask = mask + 1; 
   label_mask = zeros(size(mask)); 
   for i = 1: layer_num
        label_mask(mask==i) = layers(i).label; 
   end
   
end

