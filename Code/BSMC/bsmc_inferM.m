function state = bsmc_inferM(state, options)

% Given the motion vectors for the background objects, infer the optical
% flow field. Use gaussian belief propagation GaBP. 

% Compute Vxx
h = size(state.frame,1);
w  = size(state.frame,2);
% n = w * h;
Vxx = bsmc_computeVxx(w,h);
Vxx = Vxx ./ (options.sig_edge^2);

layer_num = length(state.layers); 
layers = state.layers; 

for i = 1: layer_num
    sparse_points_m = layers(i).traj_store.sparse_points_m; 
    
    if(layers(i).label == state.bgClust) % if it's background. 
        % Background is pretty rigidm no need to compute var
        [M_b, sig2_b] = computeFlow(Vxx, [h w], sparse_points_m, options, false);
        layers(i).motion_model = struct('mean', M_b, 'var',sig2_b);
    
    else
        if (isempty(sparse_points_m))
            M_f = zeros(h,w,2);
            sig2_f = 3^2 * ones(h,w);
        else
            [M_f,sig2_f] = computeFlow(Vxx, [h w], sparse_points_m, options, false); % original true
        end
        layers(i).motion_model = struct('mean', M_f, 'var',sig2_f);
    end
end
state.layers = layers;

end

function [flow, sig2] = computeFlow(Vxx_in, sz, points, options, compute_sig2)
    %%% compute_sig2: true,  compute var otherwise not 
    
    is = sub2ind(sz, points(2,:), points(1,:));
    n = prod(sz);
    m = size(points,2);
    h = sz(1);
    w = sz(2);

    % Update Vxx
    Vxx = Vxx_in;
    VxxDiff = sparse(is,is,repmat((1/(options.sig_motion^2)), length(is),1),n,n);
    %Vxx(sub2ind([n n], is,is)) = Vxx(sub2ind([n n], is, is)) + (1/(options.sig_motion^2));
    Vxx = Vxx + VxxDiff;

    % Compute Vxy 
    Vxy = sparse(n,m);
    Vxy(sub2ind([n m], is, 1:m)) = -(1/(options.sig_motion^2));

    % Now add labeled information.
    % Find motion vectors
    y1 = (points(3,:) - points(1,:))';
    y2 = (points(4,:) - points(2,:))';

    flow1 = Vxx\(-Vxy * y1);
    flow2 = Vxx\(-Vxy * y2);
    flow(:,:,1) = reshape(full(flow1),h,w);
    flow(:,:,2) = reshape(full(flow2),h,w);         

    if (~compute_sig2)
        sig2 = ones(h,w);
    else
        % Down scale input 8 times
        h = sz(1);
        w  = sz(2);
        if (mod(h,8) == 0 && mod(w,8) == 0)
            r = 8;
        else
            if (mod(h,5) == 0 && mod(h,5) == 0)
                r = 5;
            else
                disp('Cannot find a suitable down-scale for the image');
                sig2 = 1*ones(h,w);
                return;
            end
        end
        sig_edge = options.sig_edge * r;
        points = ceil(points/r);
        w = w/r; h = h/r;
        sz = [h w];

        % Compute sig2
        Vxx = bsmc_computeVxx(w,h);
        Vxx = Vxx ./ (sig_edge^2);
        is = sub2ind(sz, points(2,:), points(1,:));
        n = prod(sz);
        VxxDiff = sparse(is,is,repmat((1/(options.sig_motion^2)), length(is),1),n,n);
        Vxx = Vxx + VxxDiff;
        invVxx = inv(Vxx);
        sig2_downsampled = reshape(full(diag(invVxx)),h,w);
        [x,y ] = meshgrid(1:r*w,1:r*h);
        x = ceil(x/r);
        y = ceil(y/r);
        sig2 = sig2_downsampled(sub2ind(size(sig2_downsampled), y(:),x(:)));
        sig2 = reshape(sig2, r*h,r*w);
    end
end