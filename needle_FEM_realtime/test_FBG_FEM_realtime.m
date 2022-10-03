% testfile
clear;
clc;

% close all pnet connection


addpath ../sm130_interrogator_matlab/
addpath ../rawdata_process/
addpath ./WYZ_FEM
close_pnet();
filename = fullfile(tempdir,'communicate.dat');
a = exist(filename);
curvature = zeros(4,2);
[ds, ks, sb, h, l, PropertyTable] = FBG_integrated_FEM_Ogden_UTru(curvature(:,2));

if a == 0
    fid = fopen(filename,'w');
    fwrite(fid,ds,'double');
    fwrite(fid,ks,'double');
    fwrite(fid,sb,'double');
    fwrite(fid,h,'double');
    fwrite(fid,l,'double');
    fwrite(fid,curvature,'double');
    fclose(fid);
end

% get output size
sz_ds = size(ds);
sz_ks = size(ks);
sz_sb = size(sb);
sz_h = size(h);
sz_l = size(l);
sz_curvature = size(curvature);

m = memmapfile(filename, 'Writable',true, 'Format', ...
    {'double', sz_ds, 'ds';
    'double', sz_ks, 'ks';
    'double', sz_sb, 'sb';
    'double', sz_h, 'h';
    'double', sz_l, 'l';
    'double', sz_curvature, 'curvature';
    });

NumChannel = 2;
NumAA = 4;
RefData = [];
RawData = [];


% run ini_interrogator.m
interrogator = ini_interrogator('IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);

% read the reference data
RefData = mean(Read_interrogator(20,2,4,interrogator));

%while 1
while 1
    tic
    
    RawData = Read_interrogator(1,NumChannel,NumAA,interrogator);
    %toc about 0.01s
    curvature = data_process(RawData,RefData,NumChannel,NumAA);
    %toc about 0.01s
    [ds, ks, sb, h, l, PropertyTable] = FBG_integrated_FEM_Ogden_UTru(curvature(:,2));
    %toc about 0.09s

    % propertytable interval mut alphaT GammaT
    m.Data.ds(sz_ds(1),1:sz_ds(2)) = ds;
    m.Data.ks(sz_ks(1),1:sz_ks(2)) = ks;
    m.Data.sb(sz_sb(1),sz_sb(2)) = sb;
    m.Data.h(sz_h(1),sz_h(2)) = h;
    m.Data.l(sz_l(1),sz_h(2)) = l;
    m.Data.curvature(1:sz_curvature(1),1:sz_curvature(2)) = curvature;
    toc % about 0.09s
end


