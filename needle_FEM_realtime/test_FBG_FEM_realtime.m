% testfile
clear;
clc;

addpath ../sm130_interrogator_matlab/
addpath ../rawdata_process/
addpath ./WYZ_FEM

filename = fullfile(tempdir,'communicate.dat');
a = exist(filename);
[ds, ks, sb, h, l, PropertyTable] = FBG_integrated_FEM_Ogden_UTru([0;0;0;0]);

if a == 0
    fid = fopen(filename,'w');
    fwrite(fid,ds,'double');
    fwrite(fid,ks,'double');
    fwrite(fid,sb,'double');
    fwrite(fid,h,'double');
    fwrite(fid,l,'double');
    fclose(fid);
end

% get output size
sz_ds = size(ds);
sz_ks = size(ks);
sz_sb = size(sb);
sz_h = size(h);
sz_l = size(l);

m = memmapfile(filename, 'Writable',true, 'Format', ...
    {'double', sz_ds, 'ds';
    'double', sz_ks, 'ks';
    'double', sz_sb, 'sb';
    'double', sz_h, 'h';
    'double', sz_l, 'l';
    });

NumChannel = 2;
NumAA = 4;
RefData = [];
RawData = [];
curvature = zeros(4,2);

% read the reference data
RefData = mean(Read_interrogator(10,2,4,'IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1));


%while 1
while 1
    %tic
    
    RawData = Read_interrogator(1,NumChannel,NumAA,'IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);
    
    curvature = data_process(RawData,RefData,NumChannel,NumAA);
    
    [ds, ks, sb, h, l, PropertyTable] = FBG_integrated_FEM_Ogden_UTru(curvature(:,1));
    % ds 1*166 double
    % ks 1*166 double
    % sb 1*1
    % h 1*1
    % l 1*1
    % propertytable interval mut alphaT GammaT
    m.Data.ds(sz_ds(1),1:sz_ds(2)) = ds;
    m.Data.ks(sz_ks(1),1:sz_ks(2)) = ks;
    m.Data.sb(sz_sb(1),sz_sb(2)) = sb;
    m.Data.h(sz_h(1),sz_h(2)) = h;
    m.Data.l(sz_l(1),sz_h(2)) = l;
    %toc % about 0.10s
end


