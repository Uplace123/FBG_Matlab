function curvature = data_process(RawData,RefData,NumChannel,NumAA)
%load Cal_matrix.mat %for three channel needle
%load Cal_matrix_2CH.mat %for two channel needle
% cal_matrix = [ C1; C2; C3; C4] each C with a dimension of 3*2

% output
curvature = zeros(4,2);

if NumChannel == 2
    load Cal_matrix_2CH.mat
elseif NumChannel == 3
    load Cal_matrix.mat
end

row_num = NumChannel; 

% get the curvature
diff_val_mean = RawData - RefData;

for i = 1:1:NumAA
    AAindex = i:NumAA:NumAA*NumChannel;
    get_curvature = diff_val_mean(1,AAindex)* cal_matrix(i*row_num - (row_num-1) : i*row_num,:);
    curvature(i,:) = get_curvature;
end
