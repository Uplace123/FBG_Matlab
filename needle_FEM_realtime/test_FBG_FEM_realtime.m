% testfile
clear;
clc;

% close all pnet connection


addpath ../sm130_interrogator_matlab/
addpath ../rawdata_process/
addpath ./WYZ_FEM

close_pnet();

% initial data
NumChannel = 2;
RefData = [];
RawData = [];


% FEM needed input
%FEM_params;
load FEM_params.mat
NumAA = 4; % this should keep tha same as needle
curvatures = zeros(NumAA,2); %num_AA * 2
curvatures_xy = curvatures(:,1);
curvatures_xz = curvatures(:,2);

[ds, ks, xs] = FBG_integrated_FEM_Ogden_UTru(curvatures_xz,AA_lcn,...
    EBC_cur,EBC,max_outer_iter,...
    max_inner_iter,EBC_converged,...
    load_ratio,nDOF,Db,Kb,d,...
    ti, E, I, PropertyTable,sb,LM,tol,...
    outer_iter,converged,freeDOF,Nel,h);

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




% run ini_interrogator.m
interrogator = ini_interrogator('IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);

% read the reference data
RefData = mean(Read_interrogator(20,2,NumAA,interrogator));

%while 1
while 1
    tic
    
    RawData = Read_interrogator(1,NumChannel,NumAA,interrogator);
    %toc about 0.01s
    curvatures = data_process(RawData,RefData,NumChannel,NumAA);
    curvatures_xy = curvatures(:,1);
    curvatures_xz = curvatures(:,2);
    %toc about 0.01s
    [ds, ks, xs] = FBG_integrated_FEM_Ogden_UTru(curvatures_xz,AA_lcn,...
    EBC_cur,EBC,max_outer_iter,...
    max_inner_iter,EBC_converged,...
    load_ratio,nDOF,Db,Kb,d,...
    ti, E, I, PropertyTable,sb,LM,tol,...
    outer_iter,converged,freeDOF,Nel,h);
    %toc about 0.09s

    % propertytable interval mut alphaT GammaT
    m.Data.ds(sz_ds(1), 1:sz_ds(2)) = ds;
    m.Data.ks(sz_ks(1), 1:sz_ks(2)) = ks;
    m.Data.xs(sz_xs(1), 1:sz_xs(2)) = xs;
    m.Data.curvature(1:sz_curvatures(1),1:sz_curvatures(2)) = curvatures;
    toc % about 0.09s
end


