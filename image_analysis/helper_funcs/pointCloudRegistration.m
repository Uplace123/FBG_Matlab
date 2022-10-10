function [R, t, sigma] = pointCloudRegistration(vec_a, vec_b)
% 3D point clouds to point clouds registration method
% Inputs: vec_a, vec_b are cell arrays of matching 3x1 coordinates of size N;
% Outputs: 3x3 rotation matrix R, a 3x1 translation vector t, and a
% magnification factor sigma such that b = sigma*R*a + t

% Yanzhou Wang
% Oct 7 2021

%% Sanity Check
if numel(vec_a) ~= numel(vec_b)
   error('Number of elements disagrees');
   R = NaN(3, 3);
   t = NaN(3, 1);
   return
else
    for i = 1:numel(vec_a)
        if all(size(vec_a{i}) == [1, 3])
            vec_a{i} = vec_a{i}';
            warning('Vector should be 3x1. Transposing vector');
        elseif all(size(vec_b{i}) == [1, 3])
            vec_b{i} = vec_b{i}';
            warning('Vector should be 3x1. Transposing vector');
        else
            error('One or more vectors are not a 3 vector');
            R = NaN(3, 3);
            t = NaN(3, 1);
            return
        end
    end
end

%% Main Method
N = numel(vec_a);

% Pre-allocating memory
vec_a_sum = zeros(3, 1);
vec_b_sum = zeros(3, 1);
vec_a_tild = cell(1, N); % Local vectors
vec_b_tild = cell(1, N);

for i = 1:N
    vec_a_sum = vec_a_sum + vec_a{i};
    vec_b_sum = vec_b_sum + vec_b{i};
end
a_bar = 1/N*vec_a_sum; % Centroid coordinates
b_bar = 1/N*vec_b_sum;

for i = 1:N
   vec_a_tild{i} = vec_a{i} - a_bar; 
   vec_b_tild{i} = vec_b{i} - b_bar;
end

sum_norm_a_i = 0;
sum_norm_b_i = 0;
for i = 1:N
   sum_norm_a_i = sum_norm_a_i + norm(vec_a_tild{i}); 
   sum_norm_b_i = sum_norm_b_i + norm(vec_b_tild{i}); 
end
sigma = sum_norm_b_i/sum_norm_a_i; % Magnification scale

% Quaternion method to solve for R
H = zeros(3,3);
H_i = zeros(3, 3);
for i = 1:N
    a_i = sigma*vec_a_tild{i};
    b_i = vec_b_tild{i};
    for l = 1:3
        for m = 1:3
            H_i(l, m) = a_i(l)*b_i(m);
        end
    end
    H = H + H_i;
end

G = zeros(4, 4);
delta = [H(2, 3) - H(3, 2); H(3, 1) - H(1, 3); H(1, 2) - H(2, 1)];
G(1, 1) = trace(H);
G(1, 2:4) = delta';
G(2:4, 1) = delta;
G(2:4, 2:4) = H + H' - trace(H)*eye(3);
[V, D] = eig(G);
[~, max_EVAL_col] = find(D == max(D(:)));
unit_quaternion = V(:, max_EVAL_col);
R = quat2rotm(unit_quaternion');
t = b_bar - sigma*R*a_bar;
end