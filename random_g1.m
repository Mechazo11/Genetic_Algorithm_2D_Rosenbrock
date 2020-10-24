function [rand_ls] = random_g1(N,min,max, n_feature)
% N-- number of candidates to genereate
% min -- min of variable range
% manx -- max of variable range
rng(0,'twister');
a = min;
b = max;
num_to_gen = N * n_feature; % Number of elements needed to generate (N by n_feature) genes for N candidates
rand_ls = ((b-a).*rand(num_to_gen,1) + a); % Column vector here
rand_ls = reshape(rand_ls,[N 2]); % Must satisfy N * n_feature = numel(rand_ls)
% i.e if 3 candiates each having 2 features then numel(rand_ls) = 6 i.e. 6
% elements
end