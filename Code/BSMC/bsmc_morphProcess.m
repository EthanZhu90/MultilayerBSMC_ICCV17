function [ state ] = bsmc_morphProcess( state, options )
%%% remove small components and fill the hole 

mask = state.lmask;
binmask = double(mask ~= state.bgClust); 
maskWithoutSmallC = double(bwareaopen(binmask, options.SmallCompThresh));

fglabel =  unique(mask(:));
fglabel = (fglabel(fglabel ~= state.bgClust))'; 

for i = fglabel
    maskWithoutSmallC((maskWithoutSmallC == 1)& (mask==i)) = i; 
end
maskWithoutSmallC(maskWithoutSmallC == 0) = state.bgClust; 

for i = fglabel
    masktoFill = maskWithoutSmallC == i; 
    maskFilled = imfill(masktoFill,'holes');
    maskWithoutSmallC(maskFilled ==1) =i; 
end

state.lmaskAfMorph = maskWithoutSmallC;
end

