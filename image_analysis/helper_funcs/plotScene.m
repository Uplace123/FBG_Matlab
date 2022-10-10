function plotScene(designGel, CTGel, designStage, CTStage, fitNeedle, CTNeedle, Rs, ts)
% Setup
figure;
hold on;
grid on;
axis equal
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title('Scene')
% view(52, 31);

% Stage origin and orientation
stage_axes_design = {[0, 0, 0]', [3, 0, 0]', [0, 3, 0]', [0, 0, 3]'};
stage_axes_vol = cell(1, numel(stage_axes_design));
for i = 1 : numel(stage_axes_design)
    stage_axes_vol{i} = Rs * stage_axes_design{i} + ts;
end
plot3(...
    [stage_axes_vol{1}(1), stage_axes_vol{2}(1)], ...
    [stage_axes_vol{1}(2), stage_axes_vol{2}(2)], ...
    [stage_axes_vol{1}(3), stage_axes_vol{2}(3)], 'r');

plot3(...
    [stage_axes_vol{1}(1), stage_axes_vol{3}(1)], ...
    [stage_axes_vol{1}(2), stage_axes_vol{3}(2)], ...
    [stage_axes_vol{1}(3), stage_axes_vol{3}(3)], 'g')

plot3(...
    [stage_axes_vol{1}(1), stage_axes_vol{4}(1)], ...
    [stage_axes_vol{1}(2), stage_axes_vol{4}(2)], ...
    [stage_axes_vol{1}(3), stage_axes_vol{4}(3)], 'b')

% Reference frame
plot3([0, 5], [0, 0], [0, 0], 'r-', 'LineWidth', 2);
plot3([0, 0], [0, 5], [0, 0], 'g-', 'LineWidth', 2);
plot3([0, 0], [0, 0], [0, 5], 'b-', 'LineWidth', 2);

% Gel Markers
plot3(designGel(:, 1), designGel(:, 2), designGel(:, 3), 'mo');
scatter3(CTGel(:, 1), CTGel(:, 2), CTGel(:, 3), 'r*');

% Stage Markers
plot3(designStage(:, 1), designStage(:, 2), designStage(:, 3), 'co');
scatter3(CTStage(:, 1), CTStage(:, 2), CTStage(:, 3), 'b*');

% Needle
scatter3(fitNeedle(:, 1), fitNeedle(:, 2), fitNeedle(:, 3), 'MarkerEdgeColor', [0.8500 0.3250 0.0980], 'Marker', '.');
scatter3(CTNeedle(:, 1), CTNeedle(:, 2), CTNeedle(:, 3), 'MarkerEdgeColor', [0.4940 0.1840 0.5560], 'Marker', 'x');
LineH(1) = plot(nan, nan, 'k.');
LegendH{1} = 'Fit';
LineH(2) = plot(nan, nan, 'kx');
LegendH{2} = 'CT';
legend(LineH, LegendH);

hold off

drawnow;
end