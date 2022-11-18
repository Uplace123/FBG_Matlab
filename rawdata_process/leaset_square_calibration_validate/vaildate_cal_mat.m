clear;

load Cal_mat_2CH_alldata.mat 
% load Cal_mat_2CH.mat % calibration matrix get from calibration.xlsx
% load Cal_mat_2CH_v2.mat 
% the difference between Cal_mat_2CH.mat and Cal_mat_2CH_v2.mat is that:
% the previous one take mean value of trials
% the Cal_mat_2CH_alldata.mat use both data in calibration.xls and validation.xls

num_CH = 2;
num_AA = 4;
num_trial_1 = 1:1:5; % for calibration.xls
num_trial_2 = 2:1:5; % for validation.xls ignore first trial

namefile = 'validation.xls'; % data used to validate the calibration matrix
num_trial = num_trial_2;
curvature = [0.25 0.8 1 1.25]; % constant curvature curve
% 
% namefile = 'calibration.xls';
% num_trial = num_trial_1;
% curvature = [0.5 1.6 2 2.5 3.2];

index = []; % record index of each AA (each row has num_CH values)


real_mat = []; % correspoding curvature for measure_mat
num_data = 0;
sum_error = zeros(1,size(curvature,2)); % record sum abs error for each data point

for i = 1:num_AA
    index = [index i:num_AA:num_CH*num_AA];
end


fig = figure('Name','validation');
set(fig, 'Position', [60, 515, 1750, 450]);
%ax = axes('Parent',fig);
subax1 = subplot(2,2,1);
subax2 = subplot(2,2,2);
subax3 = subplot(2,2,3);
subax4 = subplot(2,2,4);

txt = num2str(curvature(1));
for i = 2 : size(curvature,2)
    txt = sprintf('%s, %s',txt, num2str(curvature(i)));
end
title_str1 = sprintf('XY plane: Desire curvature: %s',txt);
title_str2 = sprintf('XZ plane: Desire curvature: %s',txt);
output_txt = sprintf('Average abs error (curvature, %s)',txt);
title(subax1,title_str1);
title(subax4,'XY plane: Desire curvature: 0');
title(subax3,title_str2);
title(subax2,'XZ plane: Desire curvature: 0');
ylabel(subax1,'curvature');
ylabel(subax2,'curvature');
ylabel(subax3,'curvature');
ylabel(subax4,'curvature');
xlabel(subax1,'count');
xlabel(subax2,'count');
xlabel(subax3,'count');
xlabel(subax4,'count');
hold(subax1,'on');
hold(subax3,'on');
hold(subax2,'on');
hold(subax4,'on');

color_type = ['r','b','k','g','y']; % color for different curvature

for i = 1:size(curvature,2) % four curves (excluding straight)
    
    curve = num2str(curvature(i));
    cl = color_type(i);
    % temp compensation
    for tri = num_trial
        % get reference reading
        sheet_name_unbent = strcat('trial',num2str(tri),'_0mm'); 
        fbg_unbent_0d = readmatrix(namefile,'Sheet',strcat(sheet_name_unbent,'_0deg'));
        fbg_unbent_90d = readmatrix(namefile,'Sheet',strcat(sheet_name_unbent,'_90deg'));
    
        sheet_name = strcat('trial',num2str(tri),'_',curve,'mm'); 
        data = readmatrix(namefile,'Sheet',strcat(sheet_name,'_0deg')) - fbg_unbent_0d;
        fbg_curve_0d = data(:,index);
        one_size = size(fbg_curve_0d, 1);
        
        predict_curvature = [ones(one_size, 1) fbg_curve_0d(:, 1:2) ones(one_size, 1) fbg_curve_0d(:, 3:4) ones(one_size, 1) fbg_curve_0d(:, 5:6) ones(one_size, 1) fbg_curve_0d(:, 7:8)] * H;

        plot(subax1,1:size(predict_curvature,1),predict_curvature(:,1),cl);
        plot(subax2,1:size(predict_curvature,1),predict_curvature(:,2),cl);
        sum_error(i) = sum_error(i) + sum([abs(predict_curvature(:,5) - curvature(i)); abs(predict_curvature(:,6) - 0)]);

        plot(subax1,1:size(predict_curvature,1),predict_curvature(:,3),cl);

        plot(subax2,1:size(predict_curvature,1),predict_curvature(:,4),cl);
        plot(subax1,1:size(predict_curvature,1),predict_curvature(:,5),cl);

        plot(subax2,1:size(predict_curvature,1),predict_curvature(:,6),cl);
        plot(subax1,1:size(predict_curvature,1),predict_curvature(:,7),cl);

        plot(subax2,1:size(predict_curvature,1),predict_curvature(:,8),cl);

        data = readmatrix(namefile,'Sheet',strcat(sheet_name,'_90deg')) - fbg_unbent_90d;
        fbg_curve_90d = data(:,index);
        predict_curvature = [ones(one_size, 1) fbg_curve_90d(:, 1:2) ones(one_size, 1) fbg_curve_90d(:, 3:4) ones(one_size, 1) fbg_curve_90d(:, 5:6) ones(one_size, 1) fbg_curve_90d(:, 7:8)] * H;
        sum_error(i) = sum_error(i) + sum([abs(predict_curvature(:,5) - 0); abs(predict_curvature(:,6) - curvature(i))]);
        
        plot(subax4,1:size(predict_curvature,1),predict_curvature(:,1),cl);
        plot(subax3,1:size(predict_curvature,1),predict_curvature(:,2),cl);
        plot(subax4,1:size(predict_curvature,1),predict_curvature(:,3),cl);
        plot(subax3,1:size(predict_curvature,1),predict_curvature(:,4),cl);
        plot(subax4,1:size(predict_curvature,1),predict_curvature(:,5),cl);
        plot(subax3,1:size(predict_curvature,1),predict_curvature(:,6),cl);
        plot(subax4,1:size(predict_curvature,1),predict_curvature(:,7),cl);
        plot(subax3,1:size(predict_curvature,1),predict_curvature(:,8),cl);


    end
end

num_data =  2* (num_trial(end)-num_trial(1) + 1)*size(predict_curvature,1);

%legend("real curvature: 0.25","real curvature: 0.8","real curvarue: 1","real curvarture: 1.25");
hold(subax1,'off');
hold(subax2,'off');
hold(subax3,'off');
hold(subax4,'off');

% print average abs error
disp(output_txt);
disp(sum_error./num_data);