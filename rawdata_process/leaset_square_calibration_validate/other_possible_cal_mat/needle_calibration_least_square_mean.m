clear;
clc;
% correct least square method for FBG calibration

num_CH = 2;
num_AA = 4;
num_trial = 5; % 0 deg
%num_trial_90deg = 5; 
namefile = "calibration.xls";
%namefile = "validation.xls"; % if use this, make sure ignore first trial

index = []; % record index of each AA (each row has num_CH values)
measure_mat = []; % trial of 0 deg, trial of 90 deg
curvature = [0.5 1.6 2 2.5 3.2]; % constant curvature curve
%curvature = [0.25 0.8 1 1.25];
real_mat = []; % correspoding curvature for measure_mat

% construct of real_mat
for i = 1:size(curvature,2)
    row_real_mat = [];
    for j = 1:num_AA
        row_real_mat = [row_real_mat [curvature(i) 0; 0 curvature(i)]];
    end
    real_mat = [real_mat; row_real_mat];
end
% disp(real_mat);
% the dim of real_mat should be 2*num_curve * (numAA*numCH)
% disp(size(real_mat)); % 10 * 8

for i = 1:num_AA
    index = [index i:num_AA:num_CH*num_AA];
end

% get measure_mat


for i = 1:size(curvature,2) % five curves (excluding straight)
    curve = num2str(curvature(i));
    trial_0d = [];
    trial_90d = [];
    % temp compensation
    for tri = 1:num_trial
        % get reference reading
        sheet_name_unbent = strcat('trial',num2str(tri),'_0mm'); 
        fbg_unbent_0d = readmatrix(namefile,'Sheet',strcat(sheet_name_unbent,'_0deg'));
        fbg_unbent_90d = readmatrix(namefile,'Sheet',strcat(sheet_name_unbent,'_90deg'));
    
        sheet_name = strcat('trial',num2str(tri),'_',curve,'mm'); 
        data = readmatrix(namefile,'Sheet',strcat(sheet_name,'_0deg')) - fbg_unbent_0d;
        fbg_curve_0d = data(:,index);
        
        trial_0d = [trial_0d ; mean(fbg_curve_0d,1)]; % dim: 1*numAA*numCH
        
        data = readmatrix(namefile,'Sheet',strcat(sheet_name,'_90deg')) - fbg_unbent_90d;
        fbg_curve_90d = data(:,index);
        trial_90d = [trial_90d ; mean(fbg_curve_90d,1)];
        
    end
    % construct measrue_mat
    % disp(trial_0d);
    measure_mat = [measure_mat; mean(trial_0d,1); mean(trial_90d,1)];
end

% disp(measure_mat);

% the dim of measure_mat should be 2*num_curve * (numAA*numCH)
% disp(size(measure_mat)); % 10 * 8 
% with format curvature1_0: AA1_ch1 AA1_ch2 AA2_ch1 ...
%             curvature1_90:
% ...

%% least square get calibration matrix H
H = pinv([ones(2*size(curvature,2),1) measure_mat]' * [ones(2*size(curvature,2),1) measure_mat]) * [ones(2*size(curvature,2),1) measure_mat]' * real_mat;
%disp(H);
%disp(size(H)); % 9 * 8

% get the error of least square
predict = [ones(size(measure_mat,1),1) measure_mat] * H;
error = real_mat - predict;
%disp(error); % looks good!
disp(mean(abs(error),1)); % average abs error for channels in each AA

% then apply the calibration matrix to validation data
save('Cal_mat_2CH.mat','H');