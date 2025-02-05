% testfile
clear;
clc;

%% Add path
addpath ../sm130_interrogator_matlab/
addpath ../needle_FEM_realtime/WYZ_FEM
addpath ./needle_FEM_realtime/WYZ_FEM/invChol/
addpath ../rawdata_process/

%% Close all pnet connection
close_pnet();
FBG_switch = 1; % turn on FBG
                % turn off FBG

%% FEM needed input
sb = 20;
l = 145;
Db = 0; % control input 1
Kb = 0; % control input 2
ti = 25;
Nel = sb + l; % 1mm elements
Mu = [0; 5e3];
Alpha = [2; 3];
Interval = {[-sb, 0]; [0, l + 10]};
needle_length = 165; % the 18G FBG needle length
%AA_lcn_base = [65, 100, 135]; % measured from base, skipping the AA near tip due to poor readings
AA_lcn_base = []; % if FBG_switch is 0, keep AA_lcn_tip as []
AA_lcn_tip = needle_length - AA_lcn_base; % measured from tip
AA_lcn = sb + l - AA_lcn_tip; % measured from sb
% save params for plotting use
save("plot_params.mat",'sb','l','Db','Kb','ti','Nel','Mu','Alpha','Interval','AA_lcn_base');

%% FBG information
NumChannel = 3;
NumAA = 4; % this should keep tha same as needle
curvatures = zeros(NumAA,2); %num_AA * 2
curvatures_xy = curvatures(:,1);
curvatures_xz = curvatures(:,2);

%% Initialization, data and memory
d_init = zeros(2*(Nel + 1), 1);
[ds, ks, xs] = FBG_integrated_FEM_Ogden_UTru(sb, l, Db, Kb, ti, Nel, Mu, ...
                                             Alpha, Interval, curvatures_xz, AA_lcn, d_init);

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
interrogator = [];
if FBG_switch == 1
    % run ini_interrogator.m
    interrogator = ini_interrogator('IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);
    % read the reference data
    RefData = mean(Read_interrogator(20,NumChannel,NumAA,interrogator));
end

%% Main loop
while 1
    
    [ds, ks, xs] = FBG_FEM_realtime(sb, l, Db, Kb, ti, Nel, Mu, Alpha, Interval,...
        NumChannel,NumAA,interrogator,RefData,AA_lcn,FBG_switch, assem_d_k(ds, ks));
    % propertytable interval mut alphaT GammaT
    
    % story data for plotting
    m.Data.ds(sz_ds(1), 1:sz_ds(2)) = ds;
    m.Data.ks(sz_ks(1), 1:sz_ks(2)) = ks;
    m.Data.xs(sz_xs(1), 1:sz_xs(2)) = xs;
    m.Data.curvature(1:sz_curvatures(1),1:sz_curvatures(2)) = curvatures;


end
