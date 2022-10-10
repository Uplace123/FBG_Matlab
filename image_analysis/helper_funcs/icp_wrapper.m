%% Wrapping ICP algorithm with random rotations to get the best result
% target and source should be Nx3 matrices

function [R, t, err_final] = icp_wrapper(target, source, err_threshold)
num_target = size(target, 1);
num_source = size(source, 1);
if num_target ~= num_source
    warning('Numbers in target and source mismatch. Performing point rejection.')
    ratio = abs(size(target, 1) - size(source, 1))/max(size(target, 1), size(source, 1));
else
    ratio = 0;
end

err = 100; % initialize error
R = eye(3); % initialize rotation
t = zeros(3, 1); % initialize translation

while err(end) > err_threshold
    if ratio ~= 0
        [R_incre, t_incre, err] = icp(target', source', 'WorstRejection', ratio);
    else
        [R_incre, t_incre, err] = icp(target', source');
    end
    
    R = R_incre*R;
    t = R_incre*t + t_incre;
    source = apply_transform(source, R_incre, t_incre);

    if err(end) > err_threshold
        % Apply some rotation to provide a different initial condition
        % aka to "unstuck" from a local minimum
        R_rand = genRandRot();
        R = R_rand*R;
        t = R_rand*t;
        source = apply_transform(source, R_rand, zeros(3, 1));
    end
end
err_final = err(end);
end

%% Helper functions
function R_rand = genRandRot()
vec3 = pi.*rand(3,1);
R_rand = Rotx(vec3(1))*Roty(vec3(2))*Rotz(vec3(3));
end

function Rx = Rotx(a)
Rx = [1, 0, 0;
    0, cos(a), -sin(a);
    0, sin(a), cos(a)];
end

function Ry = Roty(b)
Ry = [cos(b), 0, sin(b);
    0, 1, 0;
    -sin(b), 0, cos(b)];
end

function Rz = Rotz(c)
Rz = [cos(c), -sin(c), 0;
    sin(c), cos(c), 0;
    0, 0, 1];
end