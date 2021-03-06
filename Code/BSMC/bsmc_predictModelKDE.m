function [state] = bsmc_predictModelKDE(prevState, state, options)
%    addpath('/home/elqursh/Projects/Research/Libraries/Belief Propagation/gabp-src/');
    
    if (~isfield(prevState,'layers'))
        % Initialize appearance model        
        %fprintf('\nInitializing appearance models\n');
        layer_num = length(state.layers);
        for i = 1: layer_num 
            app_x = zeros(size(state.obs,1), size(state.obs,2), options.kde_n * size(state.obs,3));
            app_counts = zeros(size(state.obs,1), size(state.obs,2)); 
            %app_labels = 0.5 .* ones(size(state.obs,1), size(state.obs,2)); %% what's this? 
           
            state.layers(i).app_model = struct('x',app_x, 'counts',app_counts); 
            
            prob_prior = (1/layer_num) * ones(size(state.obs,1), size(state.obs,2)); %% need to change 
            state.layers(i).prob_model = struct('prior', prob_prior, 'post', -1 ); % post is for later
        end
        return;
    end
    
    % Predict change to appearance model from previous frame to current
    % frame.
    
    %idx = find([traj_store.label] == bgClust); % if it's fg.
    
    layer_num = length(state.layers);
    
    for i = 1: layer_num 
        A = computeA(state.layers(i).motion_model, options); % as far,for me, A is just a shift matrix.
        mlabel = state.layers(i).label;   
        idx = find([prevState.layers.label] == mlabel); % find the corresponding app model 
        %%% if you can't find the app of this new cluster, we create a new element in app struct.
        if(isempty(idx)) 
            app_x = zeros(size(state.obs,1), size(state.obs,2), options.kde_n * size(state.obs,3));
            app_counts = zeros(size(state.obs,1), size(state.obs,2)); 
            %app_labels = 0.5 .* ones(size(state.obs,1), size(state.obs,2)); %% what's this? 
           
            prev_app = struct('x',app_x, 'counts',app_counts); 
            state.layers(i).app_model = predictState(A, prev_app, options);
        else
            state.layers(i).app_model = predictState(A, prevState.layers(idx).app_model, options);
        end
    end 
    
    %%% generate the prob_models
    %%%  it seems like prob_prior is not used later
    if (options.label_prior ==1)
        error('do not set label_prior to 1');
        fgprior = predictFGSig(A2, 1-prevState.model.bgpost, state.motion_models{2}.var, options);
    else
        for i = 1: layer_num 
            mlabel = state.layers(i).label;   
            idx = find([prevState.layers.label] == mlabel); % find the corresponding app model 
           
            if(isempty(idx)) % new cluster 
                prob_prior = (1/layer_num) * ones(size(state.obs,1), size(state.obs,2));
                state.layers(i).prob_model = struct('prior', prob_prior, 'post',-1 ); % post is for later
            else
                prob_prior = predictFG(A, prevState.layers(idx).prob_model.post); % eq(10) fg prob shift to current frame.
                state.layers(i).prob_model = struct('prior', prob_prior, 'post',-1 ); % post is for later
            end
        end
    end

end

function fgprediction = predictFG(A, prev_fgpost)
    [m,n] = size(prev_fgpost);
    % Apply the mean motion (encapsulated in A) to the prev_fgpost, such
    % that we now get 0 mean distributions.
    fgprediction = reshape(A * reshape(prev_fgpost, m*n,1),m,n);
end

function proball = predictFGSig(A, prev_fgpost, sig2, options)
    [m n] = size(prev_fgpost);
    fgprediction = reshape(A * reshape(prev_fgpost, m*n,1),m,n);
    proball = bsmc_adapGridKDE(fgprediction, sig2, options.win_size2);
end


function app = predictState(A, prev_bgapp, options)

    [m n c] = size(prev_bgapp.x);
    app_x = zeros(m,n,c);
    for i=1:c
        app_x(:,:,i) = reshape(A * reshape(prev_bgapp.x(:,:,i),m*n,1),m,n);
    end
    app_counts = reshape(A * reshape(prev_bgapp.counts,m*n,1),m,n);
    app = struct('x',app_x, 'counts',app_counts);
end

function A = computeA(motionModel, options)
    % Here we assume motion model is accurate.
    [m n ~] = size(motionModel.mean);
    [x,y] = meshgrid((1:n),(1:m));
    if (~isfield(options,'motion_subpixel') || options.motion_subpixel == 0)
        shifts_x = round(motionModel.mean(:,:,1)) + x;
        shifts_y = round(motionModel.mean(:,:,2)) + y;
        outofbounds = shifts_x < 1 | shifts_x > n | shifts_y < 1 | shifts_y > m;
        shifts_x = max(min(shifts_x, n),1);
        shifts_y = max(min(shifts_y, m),1);
        vals = ones(m,n);
        if (options.borderInitNeighbor == 0)
            vals(outofbounds) = 0;
        end
        idx = sub2ind([m n], shifts_y(:), shifts_x(:));
        is = (1:m*n);
        A = sparse(is, idx, vals(:),m*n,m*n);
    else
        % Each pixel comes from 4 other pixels 
        xx = motionModel.mean(:,:,1); xx = xx(:);
        yy = motionModel.mean(:,:,2); yy = yy(:);
        xx_floor = floor(xx); xx_ceil = ceil(xx);
        yy_floor = floor(yy); yy_ceil = ceil(yy);
        
        vec_ff = [ xx - xx_floor  yy - yy_floor ]';
        vec_fc = [ xx - xx_floor  yy_ceil - yy  ]';
        vec_cf = [ xx_ceil - xx   yy - yy_floor ]';
        vec_cc = [ xx_ceil - xx   yy_ceil - yy  ]';
        
        dist_ff = sqrt(sum(vec_ff.^2));
        dist_fc = sqrt(sum(vec_fc.^2));
        dist_cf = sqrt(sum(vec_cf.^2));
        dist_cc = sqrt(sum(vec_cc.^2));
        
        sum_dist = dist_ff + dist_fc + dist_cf + dist_cc;
        weight_ff = dist_ff ./ sum_dist;
        weight_fc = dist_fc ./ sum_dist;
        weight_cf = dist_cf ./ sum_dist;
        weight_cc = dist_cc ./ sum_dist;

        shifts_x_f = max(min(xx_floor + x(:), n),1);
        shifts_y_f = max(min(yy_floor + y(:), m),1);
        shifts_x_c = max(min(xx_ceil + x(:), n),1);
        shifts_y_c = max(min(yy_ceil + y(:), m),1);
        
        same_x = shifts_x_f == shifts_x_c;
        same_y = shifts_y_f == shifts_y_c;
        weight_ff(same_x) = weight_ff(same_x) + weight_cf(same_x); weight_cf(same_x) = 0;
        weight_fc(same_x) = weight_fc(same_x) + weight_cc(same_x); weight_cc(same_x) = 0;
        
        weight_ff(same_y) = weight_ff(same_y) + weight_fc(same_y); weight_fc(same_y) = 0;
        weight_cf(same_y) = weight_cf(same_y) + weight_cc(same_y); weight_cc(same_y) = 0;
        
        js_ff = sub2ind([m n], shifts_y_f, shifts_x_f);
        js_fc = sub2ind([m n], shifts_y_c, shifts_x_f);
        js_cf = sub2ind([m n], shifts_y_f, shifts_x_c);
        js_cc = sub2ind([m n], shifts_y_c, shifts_x_c);

        is = (1:m*n);
        A = sparse(repmat(is,1,4), ...
            [js_ff; js_fc; js_cf; js_cc]', ...
            [weight_ff weight_fc weight_cf weight_cc],m*n,m*n);
        assert(all((sum(A,2) -1) <= 20*eps));
    end
end



