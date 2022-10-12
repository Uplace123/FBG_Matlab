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


%% Close all pnet connection
close_pnet();

%% FEM needed input
sb = 20;
%l = 145;
l = 45;
ti = 25;
Nel = sb + l; % 1mm elements
Mu = [0, 0];
Alpha = [2, 3];
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
xd = [1; -0.1]; % desired tip position and orientation
Kp = diag([2, 2]); % proportional gain -> too large can result in non-converging FEM
x0 = zeros(2, 1); % initial tip position and orientation
u0 = zeros(2, 1); % initial base position and orientation
ic = [x0; u0]; % initial conditions
dt = 0.1;
Db = 0; % control input 1
Kb = 0; % control input 2

% save params for plotting use
save("plot_params.mat",'sb','l','ti','Nel','Mu','Alpha','Interval','AA_lcn_base','Db','Kb');

%% Initialization, data and memory
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
% run ini_interrogator.m
interrogator = ini_interrogator('IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);
% read the reference data
RefData = mean(Read_interrogator(20,2,NumAA,interrogator));

%% Main loop
while(1)
    tic
    % get du for next step
    du = numerical_jacobian_pos_ori_control(xd,Kp,ic,sb, l, ti, Nel, Mu, Alpha,...
                                                Interval,NumChannel,NumAA,interrogator,...
                                                RefData,AA_lcn);
    % motor gain
    u = du * dt;

    disp(u(1));
    % get current needle shape, tip location
    Db = Db + u(1);
    Kb = Kb + u(2);
    [ds, ks, xs] = FBG_FEM_realtime(sb, l, Db, Kb, ti, Nel, Mu, Alpha, Interval,...
                                    NumChannel,NumAA,interrogator,RefData,AA_lcn);
    
    % update current base and tip location and orientation
    xt = [ds(end);ks(end)]; % current tip position and orientation
    ut = [ds(1);ks(1)]; % current base position and orientation
    ic = [xt;ut];

    % story data for plotting
    m.Data.ds(sz_ds(1), 1:sz_ds(2)) = ds;
    m.Data.ks(sz_ks(1), 1:sz_ks(2)) = ks;
    m.Data.xs(sz_xs(1), 1:sz_xs(2)) = xs;
    m.Data.curvature(1:sz_curvatures(1),1:sz_curvatures(2)) = curvatures;
    toc
end
