% run_test_gradient.m
%   This just calls test gradient on the dynamics of the system so that I
%   don't have to call it from the command line. 

function [grad, num_grad] = run_test_gradient(params, s0)

if nargin < 2
    s0 = state_encode0(params);
end
    
[grad, num_grad] = test_gradient( @(s, params)test_constraints(s, params), @(s,params)test_constraints_grad(s, params), s0, params);

end