clear;
clc;
addpath ../

% Calibration Matrices Generation
CH = 2; % number of channels
AA = 4; % number of channels
col = AA*CH; % number of FBG columns
C11 = 0; 
C22 =0; 
C33 = 0; 
C44 = 0;
AA1_select = 1:4:col; 
AA2_select = 2:4:col; 
AA3_select = 3:4:col; 
AA4_select = 4:4:col;

fbg_0d = zeros(5,col);
fbg_90d = zeros(5,col);

for tri = 1:5 % five trials
    
    sheet_name_unbent = strcat('trial',num2str(tri),'0mm'); 
    fbg_unbent_0d = readmatrix('calibration_yera_2CH.xls','Sheet',strcat(sheet_name_unbent,'0deg'));
    fbg_unbent_90d = readmatrix('calibration_yera_2CH.xls','Sheet',strcat(sheet_name_unbent,'90deg'));

    for i = 0:4 % five curves (excluding straight)
        switch i
            case 0
                curve='0.25';
            case 1
                curve='0.8';
            case 2
                curve='1';
            case 3
                curve='1.25';
            case 4
                curve='3.125';
        end
        sheet_name = strcat('trial',num2str(tri),curve,'mm'); 
        fbg_curve_0d = readmatrix('calibration_yera_2CH.xls','Sheet',strcat(sheet_name,'0deg')) - fbg_unbent_0d;        
        fbg_0d(i+1,:) = mean(fbg_curve_0d);
    
        fbg_curve_90d = readmatrix('calibration_yera_2CH.xls','Sheet',strcat(sheet_name,'90deg')) - fbg_unbent_90d;
        fbg_90d(i+1,:) = mean(fbg_curve_90d);
    end

    % fbg data for all active areas
    AA1_fbg = [fbg_0d(:,AA1_select);fbg_90d(:,AA1_select)];
    AA2_fbg = [fbg_0d(:,AA2_select);fbg_90d(:,AA2_select)];
    AA3_fbg = [fbg_0d(:,AA3_select);fbg_90d(:,AA3_select)];
    AA4_fbg = [fbg_0d(:,AA4_select);fbg_90d(:,AA4_select)];
    
    % curvature array
    %k = zeros(1200,2);
    v = [0.25,0.8,1.0,1.25,3.125];
    k = [v,zeros(1,5);zeros(1,5),v]';
    
    % CALIBRATION MATRIX!!
    C1 = AA1_fbg\k;
    C2 = AA2_fbg\k;
    C3 = AA3_fbg\k;
    C4 = AA4_fbg\k;

    C11 = C11 + C1;
    C22 = C22 + C2;
    C33 = C33 + C3;
    C44 = C44 + C4;

end

C1 = C11/tri
C2 = C22/tri
C3 = C33/tri
C4 = C44/tri

cal_matrix = [C1;C2;C3;C4];
save('Cal_matrix_2CH.mat',"cal_matrix");