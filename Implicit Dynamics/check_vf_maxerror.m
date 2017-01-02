function [vf_maxerror, ic_max] = check_vf_maxerror(params, k)
% Like check_vf.m except it runs many iterations and just returns the worst
% value

    n = params.n;
    m = params.m;

    if nargin == 1
        k = 10;
    end

    for j = 1:k
        x0(:,j) = (100-90).*rand(n,1) + 90;
        u0(:,j) = (10-0).*rand(m,1) + 0;
        xdot0(:,j) = (10-(0)).*rand(n,1) + 0;
        
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