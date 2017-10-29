function rgs = bsmc_RGB2rgs(rgb)
    drgb = double(rgb);
    rgs = zeros(size(drgb));
    rgs(:,:,1) = (drgb(:,:,1) + 10) .* (255 ./ (sum(drgb,3) + 30)); % r
    rgs(:,:,2) = (drgb(:,:,2) + 10) .* (255 ./ (sum(drgb,3) + 30)); % g
    rgs(:,:,3) = sum(drgb,3) ./ 3;                % s 
end
