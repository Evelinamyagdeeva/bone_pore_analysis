clc;
clear;
close all;

load('bone_cSCI_1.mat');
boneMask = logical(VonlyBone);

[H, W, Z] = size(boneMask);
hullVolume = false(H, W, Z);

[Xgrid, Ygrid] = meshgrid(1:W, 1:H);

alphaValue = 15;   

for k = 1:Z

    slice2D = boneMask(:,:,k);

    if any(slice2D(:))

        [yPts, xPts] = find(slice2D);
        points = [xPts, yPts];

        if size(poin ts,1) > 3

            shp = alphaShape(points, alphaValue);

            alpha_mask = inShape(shp, Xgrid, Ygrid);

            alpha_mask = imfill(alpha_mask,'holes');

            hullVolume(:,:,k) = alpha_mask;

        end
    end
end
volshow(boneMask);
volshow(hullVolume);

differenceVolume = hullVolume & ~boneMask; 
volshow(differenceVolume);


r = 20;
se = strel('disk', r);

result = false(size(differenceVolume));

for i = 1:size(differenceVolume, 3)
    
    slice = differenceVolume(:,:,i);
    
 
    slice = imopen(slice, se);
    
    CC = bwconncomp(slice);
    
    if CC.NumObjects > 0
        
        numPixels = cellfun(@numel, CC.PixelIdxList);
        [~, idx] = max(numPixels);
        
        cleanSlice = false(size(slice));
        cleanSlice(CC.PixelIdxList{idx}) = true;
        
        result(:,:,i) = cleanSlice;
    end
end

volshow(result);

por = hullVolume & ~boneMask; 
por = por & ~result;
volshow(por);

save('pore_cSCI_1', 'por');


