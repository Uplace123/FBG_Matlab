function dis = extract_dist_from_d(d_to_match)
num_d = length(d_to_match)/2 - 1;
dis = zeros(1, num_d + 1);
    for i = 1:(num_d + 1)
        dis(i) = d_to_match(2*i - 1);
    end
end