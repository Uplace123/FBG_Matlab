%% Apply rigid transformation to set of data points
% coord: Nx3 coordinates
% R: 3x3 rotation matrix
% t: 3x1 translation vector
% tf_coord: transformed Nx3 coordinates

function tf_coord = apply_transform(coord, R, t)
num_coord = size(coord, 1);
tf_coord = zeros(size(coord));
for i = 1:num_coord
    tf_coord(i, :) = (R*coord(i, :)' + t)';
end
end