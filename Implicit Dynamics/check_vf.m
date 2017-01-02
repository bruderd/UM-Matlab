function [grad_dfdx_error, grad_dfdu_error, grad_dfdxdot_error] = check_vf(params, x0, u0, xdot0)
% k: number of times to run test_gradient (iterations)
% params: struct containing key information about system.
%
% This must be in same folder with the following functions:
%   test_gradient, test_f_dfdx, test_f_dfdu, test_f_dfdxdot, test_dfdx,
%   test_dfdu, test_dfdxdot, vf

    n = params.n;
    m = params.m;
    
    if nargin == 1
        x0 = rand(1,n);
        u0 = rand(1,m);
        xdot0 = rand(1,n);
    end

    [ grad_dfdx, numgrad_dfdx ] = test_gradient( @(arg, params)test_f_dfdx(arg, params, u0, xdot0), @(arg, params)test_dfdx(arg, params, u0, xdot0), x0, params );
    grad_dfdx_error = max( abs( grad_dfdx(:) - numgrad_dfdx(:) ) );
    
    [ grad_dfdu, numgrad_dfdu ] = test_gradient( @(arg, params)test_f_dfdu(arg, params, x0, xdot0), @(arg, params)test_dfdu(arg, params, x0, xdot0), u0, params );
    grad_dfdu_error = max( abs( grad_dfdu(:) - numgrad_dfdu(:) ) );

    [ grad_dfdxdot, numgrad_dfdxdot ] = test_gradient( @(arg, params)test_f_dfdxdot(arg, params, x0, u0), @(arg, params)test_dfdxdot(arg, params, x0, u0), xdot0, params );
    grad_dfdxdot_error = max( abs( grad_dfdxdot(:) - numgrad_dfdxdot(:) ) );

end