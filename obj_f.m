function [fval] = obj_f(X)
% Objective function
% Change this function depending on your optimization problem

% Input is a N * n_feature vector where n_feature is number of design
% variables and N is the number of design candidates per generation

% Unpack and form row vectors for x1 and x2
x1 = X(:,1);
x2 = X(:,2);
    
% Objective function
fval = 100 .* (x2 - x1.^2).^2 + (1 - x1).^2;
end