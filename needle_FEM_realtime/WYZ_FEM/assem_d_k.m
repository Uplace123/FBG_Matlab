function dk = assem_d_k(d, k)
% Re-assemble the DOF vector from d and k vectors
dim = numel(d);
dk = zeros(2*dim, 1);
for i = 1:dim
    dk(2*i - 1:2*i) = [d(i); k(i)];
end