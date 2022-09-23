function curvature = data_process(RawData,RefData,NumChannel,NumAA)
load Cal_matrix.mat
% cal_matrix = [ C1; C2; C3; C4] each C with a dimension of 3*2

% output
curvature = [];
row_num = 3; % each cal matrix for one AA is 3 by 2.
% get the curvature
diff_val = RawData - RefData;
% temperature compensation??
diff_val_mean = mean(diff_val); %mean value of each colume

for i = 1:1:NumAA
    AAindex = i:NumAA:NumAA*NumChannel;
    get_curvature = diff_val_mean(:,AAindex)* cal_matrix(i*row_num - (row_num-1) : i*row_num,:);
    curvature = [curvature;get_curvature];
end
