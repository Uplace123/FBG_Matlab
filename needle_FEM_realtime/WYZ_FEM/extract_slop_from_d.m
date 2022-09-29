function k = extract_slop_from_d(d_to_match)
num_k = length(d_to_match)/2 - 1;
k = zeros(1, num_k + 1);
    for i = 1:(num_k + 1)
        k(i) = d_to_match(2*i);
    end
end