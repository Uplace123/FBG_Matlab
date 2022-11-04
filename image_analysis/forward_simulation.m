addpath(genpath('./image_analysis/'));
% addpath(genpath('./invChol/'))
close all

%% Load Detected
% List of files to run
list = ["~/Desktop/DICOM/20221028/00284601"
"~/Desktop/DICOM/20221028/03063091"
"~/Desktop/DICOM/20221028/06145411"
"~/Desktop/DICOM/20221028/52172151"
"~/Desktop/DICOM/20221028/54395981"
"~/Desktop/DICOM/20221028/55585561"
"~/Desktop/DICOM/20221028/59067061"
];

baseline_number = 5;
base_dataset = list(baseline_number);
subs_dataset = setdiff(list(:), base_dataset, 'stable');

% get initial condition
threshold = 15000;
[Nel, ts, baseline_y_fit] = analyzeBaselineImage(base_dataset, threshold);

% pre-allocation
sbs = zeros(size(subs_dataset));
ls = zeros(size(subs_dataset));
dbs = zeros(size(subs_dataset));
kbs = zeros(size(subs_dataset));
tis = zeros(size(subs_dataset));

ds = zeros(size(subs_dataset, 2), Nel + 1);
ds_sim = zeros(size(ds));

rcms = zeros(size(subs_dataset)); % virtual RCM below the surface
coms = zeros(size(subs_dataset)); % needle COM
coms_sim = zeros(size(subs_dataset));

tip_errs = zeros(size(subs_dataset)); % mean tip errors
shape_errs = zeros(size(subs_dataset)); % mean shape erros

% get subsequent conditions
for i = 1:length(subs_dataset)
    [sbs(i), ls(i), dbs(i), kbs(i), tis(i), ds(i, :)] = ...
        analyzeSubsequentImage(subs_dataset(i), threshold, Nel, ts, baseline_y_fit);
end

%% Simulation
% Specify which function and parameters to use
% fem = @ogden_unconstrained_eng; % using unconstrained formulation and engineering stress
fem = @ogden_unconstrained_true; % using unconstrained formulation and true stress
% fem = @ogden_constrained_true; % using constrained formulation and true stress

% parameter values from compression tests
% Alpha = [-1; -1];
% Mu = [0; 1.2715e+04];

% parameter values from paper for porcine skeletal muscles
Alpha =[-1; 8.74];
Mu = [0; 3.63e+03];

% timed
ts = zeros(1, length(subs_dataset));

for i = 1:length(subs_dataset)
    Interval = {[-sbs(i), 0]; [0; ls(i)]};
    tic
    [ds_sim(i, :), ~] = fem(sbs(i), ls(i), dbs(i), kbs(i), tis(i), Nel, Mu, Alpha, Interval); % simulate shape
    ts(i) = toc;
    tip_errs(i) = abs(ds_sim(i, end) - ds(i, end)); % tip error
    shape_errs(i) = norm(ds_sim(i, :) - ds(i, :))/Nel; % shape error
end

% Printing
for i = 1:length(subs_dataset)
    fprintf('Error @[db, kb] = [%.2f, %.2f] [shape, tip]: [%.2f, %.2f]\n', ...
        dbs(i), kbs(i), shape_errs(i), tip_errs(i));
end

%% Plotting
f = figure;
set(f, 'Position', [130, 600, 1300, 400])
hold on
cmap = colormap(hsv(length(subs_dataset)));
for i = 1:length(subs_dataset)
    t = linspace(-sbs(i), ls(i), Nel + 1);
    color = cmap(i, :);
    plot(t, ds(i, :), '.', 'Color', color, 'MarkerSize', 3); % detected shape
    plot(t, ds_sim(i, :), '-', 'Color', [color, 0.3], 'LineWidth', 1.5); % simulated shape
end

% plot(median(rcms), 0, 'k^', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); % RCM
% plot(median(coms), 0, 'kp', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); % detected COM
% plot(median(coms_sim), 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); % simulated COM

% trick for legend
PlotH(1) = plot(nan, nan, 'k:');
LegendH{1} = 'Detected Shape';
PlotH(2) = plot(nan, nan, 'k-');
LegendH{2} = 'Simulated Shape';
% PlotH(3) = plot(nan, nan, 'k^', 'MarkerFaceColor', 'k');
% LegendH{3} = 'Needle RCM';
% PlotH(4) = plot(nan, nan, 'kp', 'MarkerFaceColor', 'k');
% LegendH{4} = 'Detected COM';
% PlotH(5) = plot(nan, nan, 'ko', 'MarkerFaceColor', 'k');
% LegendH{5} = 'Simulated COM';

% legend(PlotH, LegendH, 'Location', 'NorthEastOutside');

hold off
grid on
axis equal tight
xlabel('Depth (mm)');
ylabel('Displacement (mm)');
