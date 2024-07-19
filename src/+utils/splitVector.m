function result = splitVector(vector)
% result = utils.splitVector(vector);
%
% Given a vector, produce all possible combinatinations of ways to split
% the elements into two vectors of equal length

n = length(vector);

if mod(n, 2) ~= 0 || n == 0
    error('Input vector must have even length and cannot be empty.');
end

half_length = n / 2;
all_combinations = nchoosek(1:n, half_length);
num_combinations = size(all_combinations, 1);

result = cell(num_combinations, 1);
for i = 1:num_combinations
    part1 = all_combinations(i, :);
    part2 = setdiff(1:n, all_combinations(i, :));
    result{i} = {part1, part2};
end

end