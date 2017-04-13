function dfdx = test_dfdx(arg, params, u_const, dxdt_const)

    n = params.n;
    m = params.m;
    
    if nargin == 2
        u_const = ones(m,1);
        dxdt_const = zeros(n,1);
    end
    
    % For testing dfdx
    x = arg;        
    u = u_const;  % u is fixed
    dxdt = dxdt_const;  % dxdt if fixed
    
    [~, dfdx, ~, ~] = vf(x, u, dxdt, params);

end