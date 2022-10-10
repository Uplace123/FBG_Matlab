clear
clc;

addpath ../needle_FEM_realtime/WYZ_FEM/ 
%% FEM Parameters
S.sb = 20;
S.l = 40;
S.ti = 25;
S.Nel = 20;
S.Mu = 3000;
S.Alpha = 8;
S.function = @FBG_integrated_FEM_Ogden_UTru;
%% Control and ODE Parameters
% Requires needle tip position and orientation as feedback

% Since we have something like dx = J du and y = x, where the actual
% base manipulation is the time integration of du, the idea behind the
% simulation is to let ODE45 handle the integration of both x and u

% For actual control implementation with constant time steps, the control
% is simply u = du*dt

S.xd = [1; -0.1]; % desired tip position and orientation
S.Kp = diag([2, 2]); % proportional gain -> too large can result in non-converging FEM
tspan = [0, 1]; % simulation time span
x0 = zeros(2, 1); % initial tip position and orientation
u0 = zeros(2, 1); % initial base position and orientation
ic = [x0; u0]; % initial conditions

%% Simulation and Plotting
[t, xu] = ode45(@(t, xu) sys_ode(t, xu, S), tspan, ic);


% % Plot outputs
% f1 = figure;
% plot(t, xu(:, 1), t, xu(:, 2));
% legend('Tip Pos', 'Tip Ori')
% xlabel('Time (s)')
% title('Needle Tip States')
% % Plot base maneuvers
% f2 = figure;
% plot(t, xu(:, 3), t, xu(:, 4));
% legend('Base Pos', 'Base Ori')
% xlabel('Time (s)')
% title('Needle Base Maneuvers')
% % Plot needle shapes
% f3 = figure;
% m = linspace(-S.sb, S.l, S.Nel + 1);
% plt = plot(m, rand(1, S.Nel + 1), 'k-', 'LineWidth', 2);
% xlabel('x (mm)');
% ylabel('y (mm)');
% title('Needle Shape')
% for i = 1:length(t)
% [ds, ~] = S.function(S.sb, S.l, xu(i, 3), xu(i, 4), S.ti, S.Nel, S.Mu, S.Alpha);
% axis equal
% grid on
% set(plt, 'YData', ds);
% drawnow;
% pause(0.2);
% end
%% Auxiliary Functions
% System ODE
function dxudt = sys_ode(t, xu, S)
x = xu(1:2);
u = xu(3:4);
du = sys_control(x, u, S); % du is the control inputs, but dependent on u, which is the previous needle base configuration
[ds, ks] = S.function(S.sb, S.l, du(1), du(2), S.ti, S.Nel, S.Mu, S.Alpha);
dx1 = ds(end);
dx2 = ks(end);
dx = [dx1; dx2];
dxudt = [dx; du];
end

% Control inputs based on feedback linearization and computed torque law
function du = sys_control(x, u, S)
jac = numeric_jacobian(@(u) input_output_fem(u, S), [u(1), u(2)]);
du = -inv(jac)*S.Kp*(x - S.xd);
end

% Numerical jacobian based on FEM simulation
function y = input_output_fem(u, S)
[ds, ks] = S.function(S.sb, S.l, u(1), u(2), S.ti, S.Nel, S.Mu, S.Alpha);
y = zeros(2, 1);
y(1) = ds(end);
y(2) = ks(end);
end
