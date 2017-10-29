function [state] = bsmc_updateModelKDE(state, options)    
    [~ , ~, c] = size(state.obs);
    
    layer_num = length(state.layers);
    
    %%% the other option is use the lmask beforeMorph 
    lmask = state.lmaskAfMorph;
    
    % Update bg appearance model for each layer
    for j = 1:layer_num
        mask_layer = zeros(size(lmask)); 
        mask_layer(lmask== state.layers(j).label) = 1; % mask for a single layer

        %mask = state.layers(j).mask;
        app = state.layers(j).app_model; 
        app.counts =  app.counts + mask_layer;% Counts are gauaranteed to be > 0 at this point. Since initail mask has all ones.
        ix = find(mask_layer);
        [im,in] = ind2sub(size(mask_layer),ix);
        idx = mod(app.counts(ix)-1, options.kde_n) + 1; % Must assume that counts > 0.
        for i=1:c
            ij = sub2ind(size(state.obs),im, in, repmat(i,size(im,1),1));
            ik = sub2ind(size(app.x), im, in, (idx-1)*c+ i);
            app.x(ik) = state.obs(ij);
        end
        state.layers(j).app_model = app; 
    end
    
    
    
    
% %     % Update bg appearance model
% %     bgapp.x = state.model.bgapp_p.x;
% %     bgapp.counts = state.model.bgapp_p.counts + state.seg.mask{1}; % Counts are gauaranteed to be > 0 at this point. Since initail mask has all ones.
% %     ix = find(state.seg.mask{1});
% %     [im,in] = ind2sub(size(state.seg.mask{1}),ix);
% %     idx = mod(bgapp.counts(ix)-1, options.kde_n) + 1; % Must assume that counts > 0.    
% %     for i=1:c
% %         ij = sub2ind(size(state.obs),im, in, repmat(i,size(im,1),1));
% %         ik = sub2ind(size(bgapp.x), im, in, (idx-1)*c+ i);
% %         bgapp.x(ik) = state.obs(ij); 
% %     end
% %     
% %     % Update foreground appearance model
% %     fgapp.x = state.model.fgapp_p.x;
% %     fgapp.counts = state.model.fgapp_p.counts + state.seg.mask{2}; % Counts are gauaranteed to be > 0 at this point. Since initail mask has all ones.
% %     ix = find(state.seg.mask{2});
% %     [im,in] = ind2sub(size(state.seg.mask{2}),ix);
% %     idx = mod(fgapp.counts(ix)-1, options.kde_n) + 1; % Must assume that counts > 0.    
% %     for i=1:c
% %         ij = sub2ind(size(state.obs),im, in, repmat(i,size(im,1),1));
% %         ik = sub2ind(size(fgapp.x), im, in, (idx-1)*c+ i);
% %         fgapp.x(ik) = state.obs(ij); 
% %     end
end
