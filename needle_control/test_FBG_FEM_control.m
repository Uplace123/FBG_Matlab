% test needle control which recieve real time FEM feedback

clear;
clc;

%% Add path
% for interrogator reading and data process
addpath ../sm130_interrogator_matlab/
addpath ../rawdata_process/
% for FEM
addpath ../needle_FEM_realtime/WYZ_FEM/invChol/
addpath ../needle_FEM_realtime/WYZ_FEM/
addpath ../needle_FEM_realtime/ 
% for needle control

% for motor control
addpath ./Galil_MATLAB_API/

%% interrogator and Galil motor controller params
close_pnet();
interrogator_ip = '192.168.1.11';
interrogator_port = 1852;
motor_controller_ip = '192.168.1.201';

%% FEM needed input
sb = 27.087;
%l = 145;
l = 50.91;
ti = 25;
Nel = round(sb + l); % 1mm elements
Mu = [0, 1.2715e4];
%Mu = [0,0];
Alpha = [2, -1];
Interval = {[-sb, 0], [0, l + 10]};
needle_length = 165; % the 18G FBG needle length
%AA_lcn_base = [65, 100, 135]; % measured from base, skipping the AA near tip due to poor readings
AA_lcn_base = [];
AA_lcn_tip = needle_length - AA_lcn_base; % measured from tip
AA_lcn = sb + l - AA_lcn_tip; % measured from sb

%% FBG information
NumChannel = 2;
NumAA = 4; % this should keep tha same as needle

%% control and ode parameters
xd = [0; 0.1]; % desired tip position and orientation
Kp = diag([2, 2]); % proportional gain -> too large can result in non-converging FEM
x0 = zeros(2, 1); % initial tip position and orientation
u0 = zeros(2, 1); % initial base position and orientation
ic = [x0; u0]; % initial conditions
dt = 0.02;
Db = 0; % control input 1
Kb = 0; % control input 2

% save params for plotting use
save("plot_params.mat",'sb','l','ti','Nel','Mu','Alpha','Interval','AA_lcn_base','Db','Kb');

%% Initialization, data and memory
% run ini_interrogator.m
interrogator = ini_interrogator('IPaddress',interrogator_ip,'Port',interrogator_port,'ReadTimeout',0.1);
% run ini_motor_controller.m
g = ini_motor_controller(motor_controller_ip);
curvatures = zeros(NumAA,2); %num_AA * 2
curvatures_xy = curvatures(:,1);
curvatures_xz = curvatures(:,2);
[ds, ks, xs] = FBG_integrated_FEM_Ogden_UTru(sb, l, Db, Kb, ti, Nel, Mu, ...
                                             Alpha, Interval, curvatures_xz, AA_lcn);

filename = fullfile(tempdir,'communicate.dat');
a = exist(filename, 'file');
if a == 0
    %plot needed parameters
    fid = fopen(filename,'w');
    fwrite(fid,ds,'double');
    fwrite(fid,ks,'double');
    fwrite(fid,xs,'double');
    fwrite(fid,curvatures,'double');
    fclose(fid);
end

% get output size
sz_ds = size(ds);
sz_ks = size(ks);
sz_xs = size(xs);
sz_curvatures = size(curvatures);

m = memmapfile(filename, 'Writable',true, 'Format', ...
    {'double', sz_ds, 'ds';
    'double', sz_ks, 'ks';
    'double', sz_xs, 'xs';
    'double', sz_curvatures, 'curvature';
    });

%% Read interrogator data for FBG initialization
RefData = [];
% read the reference data
RefData = mean(Read_interrogator(20,2,NumAA,interrogator));
total_y = 0;
total_k = 0;
total_x = 0;
%% Main loop
while(1)
    %tic
    % get du for next step
    du = numerical_jacobian_pos_ori_control(xd,Kp,ic,sb, l, ti, Nel, Mu, Alpha,...
                                            Interval,NumChannel,NumAA,interrogator,...
                                             RefData,AA_lcn);
    % motor gain
    [dx,dy,dr,Db,Kb] = robot_geometric(du(1)*dt,du(2)*dt,Db,Kb);
    total_y = total_y + dy;
    total_k = total_k + dr;
    total_x = total_x + dx;
    % control motor
    % current dont use dk
    Input_RelPos_X = -round(dx*1000);
    Input_RelPos_Y = -round(dy*1000);
    Input_Rotation = round(dr*7031.25);
    give_pos=strcat('PR ,,',num2str(Input_RelPos_Y),',',  num2str(Input_RelPos_X), ',',num2str(Input_Rotation));
    galil_command(g, give_pos);
    galil_command(g, 'BG CDE');
    galil_command(g, 'MC CDE');
    
    [ds, ks, xs] = FBG_FEM_realtime(sb, l, Db, Kb, ti, Nel, Mu, Alpha, Interval,...
                                    NumChannel,NumAA,interrogator,RefData,AA_lcn);
    
    % update current base and tip location and orientation
    xt = [ds(end);ks(end)]; % current tip position and orientation
    ut = [ds(1);ks(1)]; % current base position and orientation
    ic = [xt;ut];
    % calculate error
    e = norm(xd - xt);
    disp("error:");
    disp(e);
    if e <= 1e-3
        disp("done!Enter to go back!");
        %disp(total_y);
        %disp(total_k);
        %disp(total_x);
%         disp(Db);
%         disp(Kb);
        pause();
        % go back to home
        galil_command(g, 'PA ,,0,0,0');
        galil_command(g, 'BG CDE');
        galil_disconnect(g);
        clear g;
        break
    end

    % story data for plotting
    m.Data.ds(sz_ds(1), 1:sz_ds(2)) = ds;
    m.Data.ks(sz_ks(1), 1:sz_ks(2)) = ks;
    m.Data.xs(sz_xs(1), 1:sz_xs(2)) = xs;
    m.Data.curvature(1:sz_curvatures(1),1:sz_curvatures(2)) = curvatures;
    %toc
end

clear g
