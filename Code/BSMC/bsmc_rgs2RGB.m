function rgb = bsmc_rgs2RGB(rgs)
    rgb = zeros(size(rgs));
    rgb(:,:,1) = (rgs(:,:,1) .* ( (3.* rgs(:,:,3) + 30) ./ 255)) - 10;
    rgb(:,:,2) = (rgs(:,:,2) .* ( (3.* rgs(:,:,3) + 30) ./ 255)) - 10;
    rgb(:,:,3) = 3 .* rgs(:,:,3) - (rgb(:,:,1) + rgb(:,:,2));
end