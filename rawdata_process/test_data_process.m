%%
clear;
clc;
addpath ../sm130_interrogator_matlab/
NumChannel = 3;
NumAA = 4;
RefData = [];
RawData = [];
curvature = zeros(4,2);
sheet_name_unbent = strcat('Curve_0_Trial_',num2str(5));
RefData = readmatrix('Calibration_FBG_Data_Collection_0d.xls','Sheet',sheet_name_unbent);

% using data bags from Alex
for tr = 1:5
    sheet_name = strcat('Curve_2_Trial_' , num2str(tr));
    RawData = readmatrix('Validation_FBG_Data_Collection_0d.xls','Sheet',sheet_name);
    curvature = curvature + data_process(RawData,RefData,NumChannel,NumAA);
end
curvature = curvature ./ tr;
disp(curvature);
 
%% using data from interrogator
while 1
    RawData = Read_interrogator(1,2,4,'IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);
    RawData(:,9:12) = RawData(:,1:4);
    curvature = data_process(RawData,RefData,NumChannel,NumAA);
    disp(curvature);
end