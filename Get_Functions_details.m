function [lb, ub, dim, fobj] = Get_Functions_details(F_name)
% GET_FUNCTIONS_DETAILS Returns details of a fixed test function
%   Input: F_name - Function name (ignored, not used)
%   Output: lb - Lower bound vector
%           ub - Upper bound vector
%           dim - Dimension
%           fobj - Objective function handle

% Use simple Sphere function as the test function
dim = 30;
lb = -100 * ones(1, dim);
ub = 100 * ones(1, dim);
fobj = @(x) sum(x.^2);
end
