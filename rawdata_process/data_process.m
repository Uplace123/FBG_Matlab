function curvature = data_process(RawData,RefData,NumChannel,NumAA)
% modified in 11/18/2022
%load  %for three channel needle
%load Cal_mat_2CH_alldata.mat %for two channel needle
% cal_matrix H = (num_AA*num_CH + 1) * 2*NumAA each C with a dimension of 3*2

% output
curvature = zeros(2,NumAA);

if NumChannel == 2
    load Cal_mat_2CH_alldata.mat
elseif NumChannel == 3
    load Cal_matrix.mat
end

row_num = NumChannel; 

% get the curvature
diff_val_mean = mean(RawData,1) - RefData;

get_curvature = [1 diff_val_mean] * H;

index = [[1:2:NumAA*2]' [2:2:NumAA*2]'];
curvature = get_curvature(index);
