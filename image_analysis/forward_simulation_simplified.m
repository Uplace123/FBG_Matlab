addpath(genpath('./helper_funcs/'));
addpath(genpath('../needle_FEM_realtime/WYZ_FEM/invChol/'))
close all

%% Load Detected
% List of files to run
list = ["~/Desktop/DICOM/20221104/00458911"
"~/Desktop/DICOM/20221104/03301001"
"~/Desktop/DICOM/20221104/05268201"
"~/Desktop/DICOM/20221104/07024091"
"~/Desktop/DICOM/20221104/08488471"
"~/Desktop/DICOM/20221104/59007441"
"~/Desktop/DICOM/20221104/70201291"
];

threshold = 15000;

% pre-allocation

xps = cell(numel(list), 1);
yps = cell(numel(list), 1);
zps = cell(numel(list), 1);

% get subsequent conditions
for i = 1:numel(list)
    [xp, yp, zp] = analyzeImage(list(i), threshold);
    xps{i} = xp;
    yps{i} = yp;
    zps{i} = zp;
end

% find the reference curve
index = 0;
base_value = 10000;
offset = 0;
for i = 1:numel(list)
    if abs(yps{i}(1)) < base_value
        index = i;
        base_value = abs(yps{i}(1));
        offset = yps{i}(1);
    end
end


% each curve 
for i = 1:numel(list)
    yps{i} = yps{i} - offset;
end



%% Plotting
f = figure;
set(f, 'Position', [130, 600, 1300, 400])
hold on
cmap = colormap(hsv(numel(list)));
for i = 1:numel(list)
    color = cmap(i, :);
    plot(xps{i}, yps{i}, '-', 'Color', color, 'MarkerSize', 3); % detected shape
end

hold off
grid on
axis equal tight
xlabel('Depth (mm)');
ylabel('Displacement (mm)');