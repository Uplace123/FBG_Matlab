clear;
clc;
% Calibration Matrices Generation
% Before Running this segment: Obtain Excel File Containing FBG Data
% This will generate the #AAs Calibration Matrices based on 
% 'Calibration_FBG_Data_Collection_0d.xls'
% 'Calibration_FBG_Data_Collection_90.xls'
addpath ../
CH = 3; % number of channels 
col = 4*CH; % number of AA, columns
AA1_select = 1:4:col; 
AA2_select = 2:4:col; 
AA3_select = 3:4:col; 
AA4_select = 4:4:col;

% File Name: Make sure this is the one you are using!
cal_0_raw = 'Calibration_FBG_Data_Collection_0d.xls'; % Calibration 0d File name!
cal_90_raw = 'Calibration_FBG_Data_Collection_90d.xls'; % Calibration 90d File name!

% initialization
C11 = 0; 
C22 =0; 
C33 = 0; 
C44 = 0;

fbg_0d = zeros(6,col);
fbg_90d = zeros(6,col);

for tri = 1:5 % five trials
    
    sheet_name_unbent = strcat('Curve_0_Trial_',num2str(tri)); 
    fbg_unbent_0d = readmatrix(cal_0_raw,'Sheet',sheet_name_unbent);
    fbg_unbent_90d = readmatrix(cal_90_raw,'Sheet',sheet_name_unbent);

    for i = 0:5 % six curves (excluding straight)
        % do temperature compensation
        sheet_name = strcat('Curve_',num2str(i+1),'_Trial_',num2str(tri));
        fbg_curve_0d = readmatrix(cal_0_raw,'Sheet',sheet_name) - fbg_unbent_0d;        
        fbg_0d(i+1,:) = mean(fbg_curve_0d);
    
        fbg_curve_90d = readmatrix(cal_90_raw,'Sheet',sheet_name) - fbg_unbent_90d;
        fbg_90d(i+1,:) = mean(fbg_curve_90d);
    end

    % fbg data for all active areas
    AA1_fbg = [fbg_0d(:,AA1_select);fbg_90d(:,AA1_select)];
    AA2_fbg = [fbg_0d(:,AA2_select);fbg_90d(:,AA2_select)];
    AA3_fbg = [fbg_0d(:,AA3_select);fbg_90d(:,AA3_select)];
    AA4_fbg = [fbg_0d(:,AA4_select);fbg_90d(:,AA4_select)];
    
    % curvature array
    %k = zeros(1200,2);
    v = [0.5,1.6,2.0,2.5,3.2,4.0];
    k = [v,zeros(1,6);zeros(1,6),v]';
    
    % CALIBRATION MATRIX!!
    % use least square method
    C1 = AA1_fbg\k;
    C2 = AA2_fbg\k;
    C3 = AA3_fbg\k;
    C4 = AA4_fbg\k;

    C11 = C11 + C1;
    C22 = C22 + C2;
    C33 = C33 + C3;
    C44 = C44 + C4;

end

% Final calibration matrix is the average over the five trials
C1 = C11/tri;
C2 = C22/tri;
C3 = C33/tri;
C4 = C44/tri;

cal_matrix = [C1;C2;C3;C4];

save('Cal_matrix.mat',"cal_matrix");

