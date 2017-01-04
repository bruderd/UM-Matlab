function [vf_maxerror, ic_max] = check_vf_maxerror(params, k)
% Like check_vf.m except it runs many iterations and just returns the worst
% value

    n = params.n;
    m = params.m;
    x_0 = params.x0;
    u_0 = params.u0;

    if nargin == 1
        k = 10;
    end

    for j = 1:k
%         x0(:,j) = (1.5-x_0).*rand(n,1) + x_0;         % simpler, sets ranges of all states. 
        x0(1,j) = (30-0).*rand(1,1) + 0;    % pressure
        x0(2,j) = (1.5-x_0(2)).*rand(1,1) + x_0(2);   % fiber angle      %only works for positive gama0
        x0(3,j) = (5*x_0(3)-x_0(3)).*rand(1,1) + x_0(3);    % radius
        x0(4,j) = (1.5*x_0(4) - 0.75*x_0(4)).*rand(1,1) + 0.75*x_0(4);  % length
        
        u0(:,j) = (100-u_0).*rand(m,1) + u_0;
        xdot0(:,j) = (10-(-10)).*rand(n,1) + -10;
        
        [grad_dfdx_error, grad_dfdu_error, grad_dfdxdot_error] = check_vf(params, x0(:,j), u0(:,j), xdot0(:,j));
        dx_error(j) = grad_dfdx_error;
        du_error(j) = grad_dfdu_error;
        dxdot_error(j) = grad_dfdxdot_error;
    end
    
    [maxdx_error, index1] = max(dx_error);
    [maxdu_error, index2] = max(du_error);
    [maxdxdot_error, index3] = max(dxdot_error);    
    
    
    vf_maxerror = [maxdx_error, maxdu_error, maxdxdot_error];
    
% This is from when there used to be 3 separate outputs. Now they are consolodated into a 1x3 vector    
%     grad_dfdx_maxerror = maxdx_error;
%     grad_dfdu_maxerror = maxdu_error;
%     grad_dfdxdot_maxerror = maxdxdot_error;
    
    ic_dx = [x0(:,index1)', u0(:,index1)', xdot0(:,index1)']';
    ic_du = [x0(:,index2)', u0(:,index2)', xdot0(:,index2)']';
    ic_dxdot = [x0(:,index3)', u0(:,index3)', xdot0(:,index3)']';
    
    ic_max = [ic_dx, ic_du, ic_dxdot];  % these are the randomly generated initial conditions for the worst cases
    
end