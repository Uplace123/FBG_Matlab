%%
clear;
clc;
addpath ../sm130_interrogator_matlab/
NumChannel = 2;
NumAA = 4;
RefData = [];
RawData = [];
curvature = zeros(4,2);
%%

% get reference data
sheet_name_unbent = strcat('trial',num2str(5),'0mm'); 
RefData = readmatrix('calibration_yera_2CH.xls','Sheet',strcat(sheet_name_unbent,'0deg'));

% using data bags from Alex
for tr = 1:5
    sheet_name = strcat('trial',num2str(tr),'1.6','mm'); 
    RawData = readmatrix('calibration_val_yera_2CH.xls','Sheet',strcat(sheet_name,'0deg'));
    curvature = curvature + data_process(RawData,RefData,NumChannel,NumAA);
end
curvature = curvature ./ 5;
disp(curvature);
 
%% using data from interrogator

% read the reference data
RefData = mean(Read_interrogator(10,2,4,'IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1));
while 1
    tic
    RawData = Read_interrogator(1,NumChannel,NumAA,'IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);
    %disp(RawData-RefData);
    %disp(RefData);
    
    toc
    curvature = data_process(RawData,RefData,NumChannel,NumAA);
    
    disp(curvature);
end