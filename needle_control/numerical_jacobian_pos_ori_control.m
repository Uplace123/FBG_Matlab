function du = numerical_jacobian_pos_ori_control(xd, Kp, ic, sb, l, ti, Nel, Mu, Alpha, Interval, NumChannel, NumAA, interrogator, RefData, AA_lcn)
%% FEM and FBG Parameters
S.sb = sb;
S.l = l;
S.ti = ti;
S.Nel = Nel;
S.Mu = Mu;
S.Alpha = Alpha;
S.function = @FBG_FEM_realtime;
S.Interval = Interval;
S.NumChannel = NumChannel;
S.NumAA = NumAA; 
S.interrogator = interrogator;
S.RefData = RefData;
S.AA_lcn = AA_lcn;

%% Control and ODE Parameters
% Requires needle tip position and orientation as feedback

% Since we have something like dx = J du and y = x, where the actual
% base manipulation is the time integration of du, the idea behind the
% simulation is to let ODE45 handle the integration of both x and u

% For actual control implementation with constant time steps, the control
% is simply u = du*dt

S.xd = xd; % desired tip position and orientation
S.Kp = Kp; % proportional gain -> too large can result in non-converging FEM
%tspan = [0, 1]; % simulation time span
% x0 = zeros(2, 1); % initial tip position and orientation
% u0 = zeros(2, 1); % initial base position and orientation
% ic = [x0; u0]; % initial conditions

%% Simulation and Plotting
%[t, xu] = ode45(@(t, xu) sys_ode(t, xu, S), tspan, ic);

%% get realtime control outputs

% output
x = ic(1:2);
u = ic(3:4);
du = sys_control(x, u, S);

%% Auxiliary Functions
% System ODE
% function du = sys_ode(xu, S)
% x = xu(1:2);
% u = xu(3:4);
% [du,ds,ks,xs] = sys_control(x, u, S); % du is the control inputs, but dependent on u, which is the previous needle base configuration
% no need to get dx
% [ds, ks] = S.function(S.sb, S.l, du(1), du(2), S.ti, S.Nel, S.Mu, S.Alpha,S.Interval,S.NumChannel,S.NumAA,S.interrogator,S.RefData,S.AA_lcn);
% dx1 = ds(end);
% dx2 = ks(end);
% dx = [dx1; dx2];
% dxudt = [dx; du];
% end

% Control inputs based on feedback linearization and computed torque law
function du = sys_control(x, u, S)
%[y,ds,ks,xs] = input_output_fem(u, S);
jac = numeric_jacobian(@(u) input_output_fem(u, S), [u(1),u(2)]);
du = -inv(jac)*S.Kp*(x - S.xd);
end

% Numerical jacobian based on FEM simulation
%function [y,ds,ks,xs] = input_output_fem(u, S)
function y = input_output_fem(u, S)
[ds, ks, ~] = S.function(S.sb, S.l, u(1), u(2), S.ti, S.Nel, S.Mu, S.Alpha,S.Interval,S.NumChannel,S.NumAA,S.interrogator,S.RefData,S.AA_lcn);
y = zeros(2, 1);
y(1) = ds(end);
y(2) = ks(end);
end

end