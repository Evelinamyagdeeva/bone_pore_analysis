clc;
clear;
close all;


matFile = 'pore_cSCI_1.mat';
vgiFile = 'cSCI_1.vgi';

if ~exist(matFile, 'file')
    error('Файл %s не найден!', matFile);
end

load(matFile); 
BW = logical(por);
clear por;

dx = 1; dy = 1; dz = 1;
if exist(vgiFile, 'file')
    txt = fileread(vgiFile);
    expr = '(?i)(voxel|resolution|sampling)[^\d\-+]*([0-9Ee\+\-\.]+)';
    tokens = regexp(txt, expr, 'tokens');
    if ~isempty(tokens)
        voxelSize = str2double(tokens{1}{2});
        dx = voxelSize; dy = voxelSize; dz = voxelSize;
        fprintf('Размер вокселя загружен из VGI: %.6f\n', dx);
    else
        warning('Параметр voxelSize не найден в VGI. Используются dx=dy=dz=1');
    end
else
    warning('VGI файл не найден. Используются физические единицы = воксели (dx=dy=dz=1)');
end

Vvox = dx * dy * dz;
AreaVox = dx * dy; 

bone = ~BW;
TV = numel(BW);
BV = sum(bone(:));
VV = sum(BW(:));

TVreal = TV * Vvox;
BVreal = BV * Vvox;
VVreal = VV * Vvox;

BVTV = BV / TV;
Porosity = VV / TV;


boneSurface = bwperim(bone, 26);
BS = sum(boneSurface(:)) * Vvox;
BSBV = BS / BVreal;

eulerNum = bweuler3d(BW, 26);
ConnD = (1 - eulerNum) / TVreal;

fprintf('\n===== ГЛОБАЛЬНЫЕ МЕТРИКИ ОБРАЗЦА =====\n');
fprintf('Общий объем (TV): %.2f\n', TVreal);
fprintf('Объем кости (BV): %.2f\n', BVreal);
fprintf('Объем пор (VV): %.2f\n', VVreal);
fprintf('Пористость (VV/TV): %.4f (%.2f%%)\n', Porosity, Porosity*100);
fprintf('BV/TV: %.4f\n', BVTV);
fprintf('Поверхность кости (BS): %.2f\n', BS);
fprintf('BS/BV: %.4f\n', BSBV);
fprintf('Число Эйлера: %d\n', eulerNum);
fprintf('Плотность связности (ConnD): %.6e\n', ConnD);


% WATERSHED
fprintf('Выполняется Watershed сегментация...\n');
D = bwdist(~BW);
D = imgaussfilt3(D, 0.5);
Dneg = -D;
Dneg(~BW) = -Inf;

mask = imextendedmax(D, 0.5);
mask = bwareaopen(mask, 10);
markers = bwlabeln(mask);

L = watershed(imimposemin(Dneg, markers));
L(~BW) = 0; 

CC = bwconncomp(L > 0, 26);
numPores = CC.NumObjects;
fprintf('Сегментация завершена. Найдено пор: %d\n', numPores);


f1 = L(1,:,:);    f2 = L(end,:,:);
f3 = L(:,1,:);    f4 = L(:,end,:);
f5 = L(:,:,1);    f6 = L(:,:,end);

borderLabels = unique([f1(:); f2(:); f3(:); f4(:); f5(:); f6(:)]);
borderLabels(borderLabels == 0) = [];


maxLabel = max(L(:));
isBorderMap = false(maxLabel, 1);
if ~isempty(borderLabels)
    isBorderMap(borderLabels) = true;
end

isOpenArray = false(numPores, 1);
for i = 1:numPores
    idx = CC.PixelIdxList{i};
    poreLabel = L(idx(1));
    if poreLabel > 0 && isBorderMap(poreLabel)
        isOpenArray(i) = true;
    end
end

fprintf('Расчет геометрии и морфометрии для каждой поры...\n');

stats = regionprops3(CC, 'Volume', 'Centroid', 'PrincipalAxisLength', 'EigenVectors', 'Extent');

idPore = (1:numPores)';
centroidX = zeros(numPores, 1);
centroidY = zeros(numPores, 1);
centroidZ = zeros(numPores, 1);
volumePore = zeros(numPores, 1);
radiusPore = zeros(numPores, 1);

extentVal = zeros(numPores, 1);
axisL1 = zeros(numPores, 1);
axisL2 = zeros(numPores, 1);
axisL3 = zeros(numPores, 1);
eigenVal1 = zeros(numPores, 1);
eigenVal2 = zeros(numPores, 1);
eigenVal3 = zeros(numPores, 1);
orientVecX = zeros(numPores, 1);
orientVecY = zeros(numPores, 1);
orientVecZ = zeros(numPores, 1);
degAnisotropy = zeros(numPores, 1);
sphericityVal = zeros(numPores, 1);
poreConnD = zeros(numPores, 1);
poreEuler = zeros(numPores, 1);
poreTypeStr = strings(numPores, 1);
angleFromZ = zeros(numPores, 1);

D_dist = bwdist(~BW);

for i = 1:numPores
    idx = CC.PixelIdxList{i};
    
   
    volReal = stats.Volume(i) * Vvox;
    volumePore(i) = volReal;
    radiusPore(i) = max(D_dist(idx)) * dx;
 
    centroidX(i) = stats.Centroid(i, 1) * dx;
    centroidY(i) = stats.Centroid(i, 2) * dy;
    centroidZ(i) = stats.Centroid(i, 3) * dz;

    extentVal(i) = stats.Extent(i);

    l1 = stats.PrincipalAxisLength(i, 1) * dx;
    l2 = stats.PrincipalAxisLength(i, 2) * dy;
    l3 = stats.PrincipalAxisLength(i, 3) * dz;
    axisL1(i) = l1; axisL2(i) = l2; axisL3(i) = l3;

    eVecs = stats.EigenVectors{i};
    if ~isempty(eVecs)
        mainDir = eVecs(:, 1);
        orientVecX(i) = mainDir(1);
        orientVecY(i) = mainDir(2);
        orientVecZ(i) = mainDir(3);
        
        cosTheta = abs(dot(mainDir, [0; 0; 1]));
        angleFromZ(i) = rad2deg(acos(min(max(cosTheta, -1), 1)));
    end
    
    eigenVal1(i) = (l1 / 2)^2;
    eigenVal2(i) = (l2 / 2)^2;
    eigenVal3(i) = (l3 / 2)^2;
  
    if l1 > 0
        degAnisotropy(i) = 1 - (l3 / l1);
    else
        degAnisotropy(i) = 0;
    end
 
    [r, c, z_idx] = ind2sub(size(BW), idx);
    rmin = max(1, min(r)-1); rmax = min(size(BW,1), max(r)+1);
    cmin = max(1, min(c)-1); cmax = min(size(BW,2), max(c)+1);
    zmin = max(1, min(z_idx)-1); zmax = min(size(BW,3), max(z_idx)+1);
    
    subPore = false(rmax-rmin+1, cmax-cmin+1, zmax-zmin+1);
    subIdx = sub2ind(size(subPore), r-rmin+1, c-cmin+1, z_idx-zmin+1);
    subPore(subIdx) = true;
    
    porePerim = bwperim(subPore, 26);
    poreAreaReal = sum(porePerim(:)) * AreaVox;
    
    if poreAreaReal > 0
        sphericityVal(i) = (pi^(1/3) * (6 * volReal)^(2/3)) / poreAreaReal;
    else
        sphericityVal(i) = 0;
    end
 
    eNum = bweuler3d(subPore, 26);
    poreEuler(i) = eNum;
    if volReal > 0
        poreConnD(i) = (1 - eNum) / volReal;
    else
        poreConnD(i) = 0;
    end
    
    if isOpenArray(i)
        poreTypeStr(i) = "Open";
    else
        poreTypeStr(i) = "Closed";
    end
end


outputTable = table(...
    idPore, ...
    centroidX, centroidY, centroidZ, ...
    radiusPore, ...
    volumePore, ...
    extentVal, ...
    axisL1, axisL2, axisL3, ...
    eigenVal1, eigenVal2, eigenVal3, ...
    orientVecX, orientVecY, orientVecZ, ...
    degAnisotropy, ...
    sphericityVal, ...
    poreEuler, ...
    poreConnD, ...
    poreTypeStr, ...
    angleFromZ, ...
    'VariableNames', { ...
        'id_pore', ...
        'Centroid_X', 'Centroid_Y', 'Centroid_Z', ...
        'Radius', ...
        'Volume_Pore', ...
        'Morphometry_Extent', ...
        'Shape_AxisL1', 'Shape_AxisL2', 'Shape_AxisL3', ...
        'Eigenvalue_1', 'Eigenvalue_2', 'Eigenvalue_3', ...
        'Orientation_X', 'Orientation_Y', 'Orientation_Z', ...
        'Degree_of_Anisotropy', ...
        'Sphericity', ...
        'Euler_Number', ...
        'Connectivity_Density', ...
        'Pore_Type', ...
        'Angle_Deviation_From_Z' ...
    });

csvFileName = 'metrics_cSCI_1.csv';
writetable(outputTable, csvFileName);
fprintf('\nГотово! Таблица с результатами сохранена в файл: %s\n', csvFileName);


function eulerNum = bweuler3d(BW3D, conn)
    % Вычисление характеристики Эйлера для 3D бинарных данных (быстрый срез)
    if nargin < 2
        conn = 26;
    end
    
    BW = padarray(logical(BW3D), [1 1 1], 0);
    

    v1 = BW(1:end-1, 1:end-1, 1:end-1);
    v2 = BW(2:end,   1:end-1, 1:end-1);
    v3 = BW(1:end-1, 2:end,   1:end-1);
    v4 = BW(2:end,   2:end,   1:end-1);
    v5 = BW(1:end-1, 1:end-1, 2:end);
    v6 = BW(2:end,   1:end-1, 2:end);
    v7 = BW(1:end-1, 2:end,   2:end);
    v8 = BW(2:end,   2:end,   2:end);

    configIndices = uint8(v1)      + uint8(v2)*2  + uint8(v3)*4  + uint8(v4)*8 ...
                  + uint8(v5)*16   + uint8(v6)*32 + uint8(v7)*64 + uint8(v8)*128 + 1;

    lut26 = [ ...
        0,  1,  1,  0,  1,  0,  0, -1,  1,  0,  0, -1,  0, -1, -1, -2, ...
        1,  0,  0, -1,  0, -1, -1, -2,  0, -1, -1, -2, -1, -2, -2, -3, ...
        1,  0,  0, -1,  0, -1, -1, -2,  0, -1, -1, -2, -1, -2, -2, -3, ...
        0, -1, -1, -2, -1, -2, -2, -3, -1, -2, -2, -3, -2, -3, -3, -4, ...
        1,  0,  0, -1,  0, -1, -1, -2,  0, -1, -1, -2, -1, -2, -2, -3, ...
        0, -1, -1, -2, -1, -2, -2, -3, -1, -2, -2, -3, -2, -3, -3, -4, ...
        0, -1, -1, -2, -1, -2, -2, -3, -1, -2, -2, -3, -2, -3, -3, -4, ...
        -1, -2, -2, -3, -2, -3, -3, -4, -2, -3, -3, -4, -3, -4, -4, -5, ...
        1,  0,  0, -1,  0, -1, -1, -2,  0, -1, -1, -2, -1, -2, -2, -3, ...
        0, -1, -1, -2, -1, -2, -2, -3, -1, -2, -2, -3, -2, -3, -3, -4, ...
        0, -1, -1, -2, -1, -2, -2, -3, -1, -2, -2, -3, -2, -3, -3, -4, ...
        -1, -2, -2, -3, -2, -3, -3, -4, -2, -3, -3, -4, -3, -4, -4, -5, ...
        0, -1, -1, -2, -1, -2, -2, -3, -1, -2, -2, -3, -2, -3, -3, -4, ...
        -1, -2, -2, -3, -2, -3, -3, -4, -2, -3, -3, -4, -3, -4, -4, -5, ...
        -1, -2, -2, -3, -2, -3, -3, -4, -2, -3, -3, -4, -3, -4, -4, -5, ...
        -2, -3, -3, -4, -3, -4, -4, -5, -3, -4, -4, -5, -4, -5, -5, -6 ...
    ] / 8;


    counts = accumarray(configIndices(:), 1, [256, 1]);
    
    if conn == 26
        eulerNum = round(sum(counts .* lut26'));
    else
        eulerNum = round(sum(counts .* (-lut26)'));
    end
end
