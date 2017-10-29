function proball = bsmc_evalKDE_padded(x, counts, obs, sig2_motion, options)
    counts = min(counts, options.kde_n);
    [m n ~] = size(obs);
    
    if (~exist('adapmed','builtin'))
        addpath('../Utils/AdaptiveMedian');
    end
    
    % Compute sig_i
    sig = zeros( size(obs));
    for c = 1:3
        channel = x(:,:,c:3:end);      
        sig(:,:,c) = reshape(adapmed(reshape(permute(abs(channel(:,:,1:(end-1)) - channel(:,:,2:end)), [3 1 2]), [], m*n), reshape(counts-1, m*n,1)), m,n);
    end
    sig = sig ./ (0.68 *sqrt(2));
 
    sig = max(sig, options.minkdesig); % 0.5 is the minimum SD to avoid a 0 SD
    sig2 = sig.^2;

    % Use appearance model to find set of labels 
    proball = zeros(m,n);
    padsize = (options.win_size -1) /2;
    padded_counts = min(padarray(counts,[padsize padsize]), options.kde_n);
    padded_x = padarray(x, [padsize padsize 0]);
    padded_sig2 = padarray(sig2, [padsize padsize 0]);

    if (options.motion_window == 1)
        tid =tic;
        for i = 0:(options.win_size -1)
            for j=0:(options.win_size - 1)
                dist2 = (padsize - i)^2 + (padsize -j)^2;
                %sig2_mask = dist2 < 4 * sig2_motion; % Sum over two standard deviations
                weights = (1./ ( 2 .* pi .* sig2_motion) ) .* exp(-dist2 ./ (2.* sig2_motion));
                prob = zeros(m, n);
                counts = padded_counts((1:m)+i,(1:n)+j,:);
                diff = zeros(m,n,3*options.kde_n);
                diff(:,:,1:3:end) = bsxfun(@minus, padded_x((1:m)+i,(1:n)+j,1:3:end),obs(:,:,1));
                diff(:,:,2:3:end) = bsxfun(@minus, padded_x((1:m)+i,(1:n)+j,2:3:end),obs(:,:,2));
                diff(:,:,3:3:end) = bsxfun(@minus, padded_x((1:m)+i,(1:n)+j,3:3:end),obs(:,:,3));
                
                diff(:,:,1:3:end) = bsxfun(@rdivide, diff(:,:,1:3:end).^2, 2 .* padded_sig2((1:m)+i,(1:n)+j,1));
                diff(:,:,2:3:end) = bsxfun(@rdivide, diff(:,:,2:3:end).^2, 2 .* padded_sig2((1:m)+i,(1:n)+j,2));
                diff(:,:,3:3:end) = bsxfun(@rdivide, diff(:,:,3:3:end).^2, 2 .* padded_sig2((1:m)+i,(1:n)+j,3));
                const = 1./ sqrt( prod( 2.* pi .* padded_sig2((1:m)+i,(1:n)+j,:) ,3));

                for k=1:options.kde_n 
                    idx = find(k <= counts);
                    p = const .* exp(-sum(diff(:,:,3*(k-1)+ (1:3)),3));
                    prob(idx) = prob(idx) + p(idx);
                end
                proball = proball + weights .* prob;
            end
        end
        fprintf(' evalkde time = %f ', toc(tid));
    else
        for i = 0:(options.win_size -1)
            for j=0:(options.win_size - 1)
                prob = zeros(m, n);
                counts = padded_counts((1:m)+i,(1:n)+j,:);
                diff = bsxfun(@minus, padded_x((1:m)+i,(1:n)+j,:),obs);
                diff2 = bsxfun(@rdivide, diff.^2, 2 .* padded_sig2((1:m)+i,(1:n)+j,:));
                const = 1./ sqrt( prod( 2.* pi .* padded_sig2((1:m)+i,(1:n)+j,:) ,3));

                for k=1:options.kde_n 
                    idx = find(k <= counts);
                    p = const .* exp(-sum(diff2(:,:,3*(k-1)+ (1:3)),3));
                    prob(idx) = prob(idx) + p(idx);
                end
                proball = max(proball, prob);
            end
        end
    end
end
