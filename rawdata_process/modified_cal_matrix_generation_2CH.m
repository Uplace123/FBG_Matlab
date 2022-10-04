clear;
clc;

% Calibration Matrices Generation modified
% use both calibration and validation data
% to generate calibration matrix in order
% to minimize the error

CH = 2; % number of channels
num_AA = 4; % number of AA
col = num_AA*CH; % number of FBG columns
C11 = 0; 
C22 = 0; 
C33 = 0; 
C44 = 0;
AA1_select = 1:num_AA:col;  
AA2_select = 2:num_AA:col; 
AA3_select = 3:num_AA:col; 
AA4_select = 4:4:col; 

fbg_0d = zeros(5,col); % 5 is trial num
fbg_90d = zeros(5,col);
filename = 'calibration.xls';
% curvature array


for tri = 1:5 % five trials
    
    sheet_name_unbent = strcat('trial',num2str(tri),'_0mm'); 
    fbg_unbent_0d = readmatrix(filename,'Sheet',strcat(sheet_name_unbent,'_0deg'));
    %fbg_unbent_0d = readmatrix(filename,'Sheet',strcat(sheet_name_unbent,'_0deg_2'));
    
    fbg_unbent_90d = readmatrix(filename,'Sheet',strcat(sheet_name_unbent,'_90deg'));
    %fbg_unbent_90d2 = readmatrix(filename,'Sheet',strcat(sheet_name_unbent,'_90deg_2'));

%     fbg_unbent_0d = [];
%     fbg_unbent_90d = [];
%     for index = 1:size(fbg_unbent_0d2,1)
%         fbg_unbent_0d(index,:) = (fbg_unbent_0d1(index,:) + fbg_unbent_0d2(index,:))/2;
%         fbg_unbent_90d(index,:) = (fbg_unbent_90d1(index,:) + fbg_unbent_90d2(index,:))/2;
%     end

    for i = 0:4 % five curves (excluding straight)
        switch i
            case 0
                curve='0.5';
            case 1
                curve='1.6';
            case 2
                curve='2';
            case 3
                curve='2.5';
            case 4
                curve='3.2';
        end
        sheet_name = strcat('trial',num2str(tri),'_',curve,'mm'); 
        fbg_curve_0d = readmatrix(filename,'Sheet',strcat(sheet_name,'_0deg')) - fbg_unbent_0d;        
        fbg_0d(i+1,:) = mean(fbg_curve_0d);    
        fbg_curve_90d = readmatrix(filename,'Sheet',strcat(sheet_name,'_90deg')) - fbg_unbent_90d;
        fbg_90d(i+1,:) = mean(fbg_curve_90d);
    end

    % fbg data for all active areas
    AA1_fbg = [fbg_0d(:,AA1_select);fbg_90d(:,AA1_select)];
    AA2_fbg = [fbg_0d(:,AA2_select);fbg_90d(:,AA2_select)];
    AA3_fbg = [fbg_0d(:,AA3_select);fbg_90d(:,AA3_select)];
    AA4_fbg = [fbg_0d(:,AA4_select);fbg_90d(:,AA4_select)];
    v = [0.5,1.6,2,2.5,3.2];
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

C1_1 = C11/tri;
C2_1 = C22/tri;
C3_1 = C33/tri;
C4_1 = C44/tri;


fbg_0d = zeros(4,col); % 4 is trial num
fbg_90d = zeros(4,col);
filename = 'validation.xls';
% curvature array

C11 = 0; 
C22 = 0; 
C33 = 0; 
C44 = 0;

for tri = 1:4 % five trials
    
    sheet_name_unbent = strcat('trial',num2str(tri),'_0mm'); 
    fbg_unbent_0d = readmatrix(filename,'Sheet',strcat(sheet_name_unbent,'_0deg'));
    
    fbg_unbent_90d = readmatrix(filename,'Sheet',strcat(sheet_name_unbent,'_90deg'));


    for i = 0:3 % four curves (excluding straight)
        switch i
             case 0
                curve = '0.25';
             case 1
                curve = '0.8';
             case 2
                curce = '1';
             case 3
                curve = '1.25';
        end
        sheet_name = strcat('trial',num2str(tri),'_',curve,'mm'); 
        fbg_curve_0d = readmatrix(filename,'Sheet',strcat(sheet_name,'_0deg')) - fbg_unbent_0d;        
        fbg_0d(i+1,:) = mean(fbg_curve_0d);    
        fbg_curve_90d = readmatrix(filename,'Sheet',strcat(sheet_name,'_90deg')) - fbg_unbent_90d;
        fbg_90d(i+1,:) = mean(fbg_curve_90d);
    end

    % fbg data for all active areas
    AA1_fbg = [fbg_0d(:,AA1_select);fbg_90d(:,AA1_select)];
    AA2_fbg = [fbg_0d(:,AA2_select);fbg_90d(:,AA2_select)];
    AA3_fbg = [fbg_0d(:,AA3_select);fbg_90d(:,AA3_select)];
    AA4_fbg = [fbg_0d(:,AA4_select);fbg_90d(:,AA4_select)];
    v = [0.25,0.8,1,1.25];
    k = [v,zeros(1,4);zeros(1,4),v]';
   
    
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

C1_2 = C11/tri;
C2_2 = C22/tri;
C3_2 = C33/tri;
C4_2 = C44/tri;

C1 = (C1_1+C1_2)/2;
C2 = (C2_1+C2_2)/2;
C3 = (C3_1+C3_2)/2;
C4 = (C4_1+C4_2)/2;

cal_matrix = [C1_1;C2_1;C3_1;C4_1]
save('Cal_matrix_2CH.mat',"cal_matrix");