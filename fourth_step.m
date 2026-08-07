clear; clc; close all;

csvFileName = 'metrics_control_1.csv';
if ~exist(csvFileName, 'file')
    error('Файл %s не найден!', csvFileName);
end

data = readtable(csvFileName);

maxPoresToPlot = 3000; 
numPores = height(data);

if numPores > maxPoresToPlot
    fprintf('Всего пор: %d. Для плавности визуализируем случайную выборку из %d пор.\n', numPores, maxPoresToPlot);
    idx = randperm(numPores, maxPoresToPlot);
    dataPlot = data(idx, :);
else
    dataPlot = data;
end

X = dataPlot.Centroid_X;
Y = dataPlot.Centroid_Y;
Z = dataPlot.Centroid_Z;

U = dataPlot.Orientation_X;
V = dataPlot.Orientation_Y;
W = dataPlot.Orientation_Z;

lengthScale = dataPlot.Shape_AxisL1 * 0.5; 
U_scaled = U .* lengthScale;
V_scaled = V .* lengthScale;
W_scaled = W .* lengthScale;


figure('Color', 'w', 'Name', ['Ориентация пор: ' csvFileName]);

hQ = quiver3(X, Y, Z, U_scaled, V_scaled, W_scaled, 0, 'LineWidth', 1.2);
hold on;

angles = dataPlot.Angle_Deviation_From_Z;
scatter3(X, Y, Z, 15, angles, 'filled');

colormap(jet); % Палитра: синий = вдоль Z, красный = перпендикулярно Z
c = colorbar;
c.Label.String = 'Отклонение от оси Z ';

xlabel('X ');
ylabel('Y ');
zlabel('Z ');
axis equal; 
grid on;
view(3);   
rotate3d on; 
