function cellmat = packageIntoCell(mat)
% Utility function that package a Nx3 matrix into N cells with 3x1 vectors

% Yanzhou Wang
% Oct 28 2021
[row, ~] = size(mat);
cellmat = cell(row, 1);
for i = 1:row
cellmat{i} = mat(i, :)';
end
end