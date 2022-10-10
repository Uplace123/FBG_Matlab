addpath(genpath('./helper_funcs/'));
addpath(genpath('../needle_FEM_realtime/WYZ_FEM/invChol/'))
close all

%% Load Detected
% FBG curvatures and locations
avg_curvatures = importfile('../write_to_file/needle_curvature_2022-10-10.txt');
XZ_curvatures = table2array(avg_curvatures(:, 5:8)); % only XZ-plane curvatures are relevant

needle_length = 165; % the 18G FBG needle length
AA_lcn_base = [65, 100, 135]; % measured from base, skipping the AA near tip due to poor readings
AA_lcn_tip = needle_length - AA_lcn_base; % measured from tip

% List of files to run
top_level_folder = "~/Desktop/DICOM/";
lower_level_folders = [
    "20221010/00382901"; 
    "20221010/02311371";
    "20221010/04109001";
    "20221010/49598581"; % baseline data
    "20221010/57094461"];

list = top_level_folder + lower_level_folders;
baseline_number = 4;
base_dataset = list(baseline_number);
subs_dataset = setdiff(list(:), base_dataset, 'stable');

% get initial condition
threshold = 9250;
[Nel, ts, baseline_y_fit] = analyzeBaselineImage(base_dataset, threshold);

% pre-allocation
sbs = zeros(size(subs_dataset));
ls = zeros(size(subs_dataset));
dbs = zeros(size(subs_dataset));
kbs = zeros(size(subs_dataset));
tis = zeros(size(subs_dataset));
intervals = cell(size(subs_dataset));
AA_lcns = zeros(numel(subs_dataset), numel(AA_lcn_base));

ds = zeros(size(subs_dataset, 2), Nel + 1);
ds_sim = zeros(size(ds));
ds_fbg_sim = zeros(size(ds));

rcms = zeros(size(subs_dataset)); % virtual RCM below the surface
coms = zeros(size(subs_dataset)); % needle COM
coms_sim = zeros(size(subs_dataset));

tip_errs = zeros(size(subs_dataset)); % mean tip errors
shape_errs = zeros(size(subs_dataset)); % mean shape erros

% get subsequent conditions
for i = 1:length(subs_dataset)
    [sbs(i), ls(i), dbs(i), kbs(i), tis(i), ds(i, :)] = ...
        analyzeSubsequentImage(subs_dataset(i), threshold, Nel, ts, baseline_y_fit);
    intervals{i} = {[-sbs(i), 0], [0, ls(i)]};
    AA_lcns(i, :) = sbs(i) + ls(i) - AA_lcn_tip; % measured from sb
end

%% Simulation
% Specify which function and parameters to use
% fem = @ogden_unconstrained_eng; % using unconstrained formulation and engineering stress
fem = @ogden_unconstrained_true; % using unconstrained formulation and true stress
% fem = @ogden_constrained_true; % using constrained formulation and true stress

fbg_fem = @FBG_integrated_FEM_Ogden_UTru; % FBG-integrated FEM model

% parameter values from compression tests
Alpha = -1;
Mu = 1.2715e+04;

% parameter values from paper for porcine
% Alpha = 8.74;
% Mu = 3.63e+03;

for i = 1:length(subs_dataset)
    [ds_sim(i, :), ~, ~] = fem(sbs(i), ls(i), dbs(i), kbs(i), tis(i), Nel, Mu, Alpha); % simulate shape
    [ds_fbg_sim(i, :), ~, ~] = fbg_fem(sbs(i), ls(i), dbs(i), kbs(i), tis(i), Nel, [0, Mu], [2, Alpha], ...
                                       intervals{i}, XZ_curvatures(i, :), AA_lcns(i, :)); % simulation with FBG
%     tip_errs(i) = abs(ds_sim(i, end) - ds(i, end)); % tip error
%     shape_errs(i) = norm(ds_sim(i, :) - ds(i, :))/Nel; % shape error
end

% % Printing
% for i = 1:length(subs_dataset)
%     fprintf('Error @[db, kb] = [%.2f, %.2f] [shape, tip]: [%.2f, %.2f]\n', ...
%         dbs(i), kbs(i), shape_errs(i), tip_errs(i));
% 
%     %% additional computations
%     % COM
%     x_spacing = (ls(i) + sbs(i))/(Nel + 1);
%     [~, com_idx] = min(abs(ds(i, :))); % detected needle COM index
%     coms(i) = com_idx*x_spacing - sbs(i);
% 
%     [~, com_sim_idx] = min(abs(ds_sim(i, :))); % simulated needle COM index
%     coms_sim(i) = com_sim_idx*x_spacing - sbs(i);
% 
%     % RCM
%     rcms(i) = -sbs(i) - dbs(i)/kbs(i);
% end
% fprintf('Mean [shape, tip] error: [%.4f, %.4f] mm\n', ...
%     mean(shape_errs), mean(tip_errs));
% fprintf('Median [detected, simulated] COM: [%.2f, %.2f] mm\n', ...
%     median(coms), median(coms_sim));
% fprintf('Median needle RCM: %.2f\n', median(rcms));

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
    plot(t, ds_fbg_sim(i, :), '--', 'Color', [color, 0.5], 'LineWidth', 1.5); % FBG + simulated shape
end

% trick for legend
PlotH(1) = plot(nan, nan, 'k:');
LegendH{1} = 'Detected Shape';
PlotH(2) = plot(nan, nan, 'k-');
LegendH{2} = 'Simulated Shape';
PlotH(3) = plot(nan, nan, 'k--');
LegendH{3} = 'FBG Simulated Shape';
legend(PlotH, LegendH, 'Location', 'northeastoutside');

hold off
grid on
axis equal tight
xlabel('Depth (mm)');
ylabel('Displacement (mm)');

% error stats
% error_stats = [shape_errs, tip_errs, coms];