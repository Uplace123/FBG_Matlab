clear;
clc;

% this script collect data from sm130 interrogator
% and process to curvature then write them into file
% in needle_rawdata, each sheet include set of rawdata read from
% interrogator. each sheet should be num_count * (num_ch*num_AA).
% in needle_curvature, include the curvature calculated from rawdata.
% each row: time curvature at AA1 AA2 AA3 AA4 as shown in header
addpath ../sm130_interrogator_matlab/
addpath ../rawdata_process/

num_ch = 2;
num_AA = 4;
num_count = 200;
filename_curvature = strcat('needle_curvature_', datestr(date,29), '.txt');
filename_rawdata = strcat('needle_rawdata_',datestr(date,29),'.xls');
header = {'time','XY_AA1','XY_AA2','XY_AA3','XY_AA4','XZ_AA1','XZ_AA2','XZ_AA3','XZ_AA4'};

close_pnet();

interrogator = ini_interrogator('IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);
% get ref data

disp("press ENTER to record the Reference data...");
pause;
sheetname = 'refdata';
RefData = mean(Read_interrogator(20,2,4,interrogator));
writematrix(RefData,filename_rawdata,'sheet',sheetname);
disp("done!");


while 1
    disp("press ENTER to record data...");
    pause;
    disp("recording..");
    sheetname = datestr(now,30);
    raw_data = Read_interrogator(num_count,num_ch,num_AA,interrogator);
    % write rawdata to xls file
    disp("write to xls file...");
    writematrix(raw_data,filename_rawdata,'sheet',sheetname);
    % get curvature num_AA * 2
    disp("write to txt file...");
    curvature = data_process(raw_data,RefData,num_ch,num_AA);
    curvature_a = {sheetname,curvature(1,1),curvature(2,1),curvature(3,1),curvature(4,1),curvature(1,2),curvature(2,2),curvature(3,2),curvature(4,2)};
    if exist(filename_curvature,'file')
        writecell(curvature_a,filename_curvature,'WriteMode','append','Delimiter','tab');
    else
        writecell(header,filename_curvature,'Delimiter','tab');
        writecell(curvature_a,filename_curvature,'WriteMode','append','Delimiter','tab');
    end
    disp("done!");

end

