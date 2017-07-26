% run_test_gradient.m
%   This just calls test gradient on the dynamics of the system so that I
%   don't have to call it from the command line. 
%
%   Uses zero as the initial condition

function [grad, num_grad] = run_test_gradient(choice, params)
  
if choice == 1
    arg0 = zeros(6,1);
elseif choice == 2
    arg0 = zeros(6,1);
elseif choice == 3
    arg0 = zeros(params.numFREEs,1);
end

[grad, num_grad] = test_gradient( @(arg, params)test_f(arg, params, choice), @(arg, params)test_fgrad(arg, params, choice), arg0, params);

end