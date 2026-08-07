clc;
clear;
close all;

vgi_file = 'control_2.vgi';
fid = fopen(vgi_file);

x=[]; y=[]; z=[];
while ~feof(fid)
    line = fgetl(fid);
    if contains(line, 'size =')
        num = sscanf(line, 'size = %d %d %d');
        x=num(1); y=num(2); z=num(3);
    end
end
fclose(fid);
clear fid line num;  

fid_vol = fopen('control_2.vol','r');
data = fread(fid_vol, x*y*z, '*uint16');
fclose(fid_vol);
clear fid_vol;      

volume = reshape(data, [x,y,z]);
volume = permute(volume, [2,1,3]);
clear data x y z;    

V = single(volume);
V = mat2gray(V);

clear volume ;

level = graythresh(V)

V = V > level;
volshow(V); 
sigma = 0.7;       
V = single(V);

edges_labeled = my_kenny(V, sigma);

Vbin = V > level;
Vbin = imfill(Vbin, 26, 'holes');
Vbin = bwareaopen(Vbin, 500);  

%volshow(Vbin);
edges = edges_labeled > 0;


se = strel('sphere', 2);
edges_fat = imdilate(edges, se);


mask = Vbin & ~edges_fat;


CC = bwconncomp(mask, 26);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx] = max(numPixels);

boneMask = false(size(mask));
boneMask(CC.PixelIdxList{idx}) = true;

%volshow(boneMask);
VonlyBone = V;
VonlyBone(~boneMask) = 0;

volshow(VonlyBone);

save('bone_cSCI_1.mat', 'VonlyBone', '-v7.3');








