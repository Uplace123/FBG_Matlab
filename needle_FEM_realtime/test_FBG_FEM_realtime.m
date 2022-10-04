% testfile
clear;
clc;

% close all pnet connection


addpath ../sm130_interrogator_matlab/
addpath ../rawdata_process/
addpath ./WYZ_FEM
close_pnet();
filename = fullfile(tempdir,'communicate.dat');
a = exist(filename, 'file');
NuNumChannel = 2;
NumAA = 4;
RefData = [];
RawData = [];

% parameters initialization
% l1
sb = 0;
% l2
l = 165;
% Constants whole length
L = sb + l;
Db = 0;
Kb = 0;
% Kb = -Db/(sb + rcm);
ti = 25;
Nel = sb + l;
Mu = 0;
Alpha = 5;

%FEM_params;
E = 250*1000; % 200GPa but in mm^2
OD = 1.27;
I = pi/4*(OD/2)^4; % in mm^4

Interval = {[-sb, 0]; [0, l]};
MuT = [0; Mu]*10^-6; % Pa, but in mm^2; 1Pa = 1e-6 N/mm^2; 1kPa = 1e-3 N/mm^2
AlphaT = [0; Alpha]; % need abs(alpha) > 1 
GammaT = [0; 0];

% PropertyTable = table(Interval, MuT, AlphaT, GammaT);
% ti = 25; % Initial length of undeformed tissue. Will affect solution

% FEM-specific constants
% Nel = 10; % Total umber of elements
Nen = 4; % Number of element DOF
DOF = 1:(2*Nel + 2);
nDOF = length(DOF);
d = zeros(nDOF, 1); % Initial guess of solution
h = L/Nel; % Finite element size
EBC = [1; 2]; % Displacement and slope of first element left node is prescribed
freeDOF = DOF;
freeDOF(EBC) = [];

% Algorithm-specific constants
max_inner_iter = 5; % Max number of iterations for Newton's method
max_outer_iter = 50; % Max number of iterations for load stepping

tol = 1e-3; % Convergence criterion

load_ratio = 1;
EBC_cur = zeros(2, 1); % For load stepping
EBC_converged = zeros(2, 1); % For load stepping from previously converged EBC

% Construction of LM
LM = zeros(Nen, Nel);
for e = 1:Nel
    LM(1, e) = 2*e - 1;
    LM(2, e) = 2*e - 0;
    LM(3, e) = 2*e + 1;
    LM(4, e) = 2*e + 2;
end

%AA_lcn = [65, 100, 135, 155];% location of AA wrt to base
AA_lcn = [65, 100, 135];
curvatures = zeros(numAA,2); %num_AA * 2
AA_crv = 1e-3 * curvature'; % curvature in XZ plane at AAs in mm
AA_er = round(AA_lcn./h) + 1; % elements where the left-moment is fixed


[ds, ks, xs, sb, h, l, PropertyTable] = FBG_integrated_FEM_Ogden_UTru(curvatures(:,2));

if a == 0
    fid = fopen(filename,'w');
    fwrite(fid,ds,'double');
    fwrite(fid,ks,'double');
    fwrite(fid,sb,'double');
    fwrite(fid,xs,'double');
    fwrite(fid,h,'double');
    fwrite(fid,l,'double');
    fwrite(fid,curvatures,'double');
    fclose(fid);
end

% get output size
sz_ds = size(ds);
sz_ks = size(ks);
sz_sb = size(sb);
sz_h = size(h);
sz_l = size(l);
sz_xs = size(xs);
sz_curvatures = size(curvatures);

m = memmapfile(filename, 'Writable',true, 'Format', ...
    {'double', sz_ds, 'ds';
    'double', sz_ks, 'ks';
    'double', sz_xs, 'xs';
    'double', sz_sb, 'sb';
    'double', sz_h, 'h';
    'double', sz_l, 'l';
    'double', sz_curvatures, 'curvature';
    });




% run ini_interrogator.m
interrogator = ini_interrogator('IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);

% read the reference data
RefData = mean(Read_interrogator(20,2,numAA,interrogator));

%while 1
while 1
    tic
    
    RawData = Read_interrogator(1,NumChannel,NumAA,interrogator);
    %toc about 0.01s
    curvatures = data_process(RawData,RefData,NumChannel,NumAA);
    %toc about 0.01s
    [ds, ks, xs, sb, h, l, PropertyTable] = FBG_integrated_FEM_Ogden_UTru(curvatures(:,2));
    %toc about 0.09s

    % propertytable interval mut alphaT GammaT
    m.Data.ds(sz_ds(1), 1:sz_ds(2)) = ds;
    m.Data.ks(sz_ks(1), 1:sz_ks(2)) = ks;
    m.Data.xs(sz_xs(1), 1:sz_xs(2)) = xs;
    m.Data.sb(sz_sb(1),sz_sb(2)) = sb;
    m.Data.h(sz_h(1),sz_h(2)) = h;
    m.Data.l(sz_l(1),sz_h(2)) = l;
    m.Data.curvature(1:sz_curvatures(1),1:sz_curvatures(2)) = curvatures;
    toc % about 0.09s
end


