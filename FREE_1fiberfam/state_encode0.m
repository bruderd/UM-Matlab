function s = state_encode0(params)
% converts states and inputs ICs into decision variable s
    x0 = params.x0;
    u0 = params.u0;
    N = params.N;
    m = params.m;
    n = params.n;

%    s = [x0; rand(N*n,1); u0; rand(N*m,1); rand(N*n,1)];  %intial decision variable s (includes x, u, and dxdt for all time steps). Initialized to all ones other than initial values

    % simpler initialization to hopefully curb nasty behavior
    for j = 1:N
        x_initial(n*(j-1) + 1 : n*(j-1) + n, 1) = x0;
    end
    
    s = [x0; x_initial; u0; zeros(N*m,1); zeros(N*n,1)];  %intial decision variable s (includes x, u, and dxdt for all time steps). Initialized to all ones other than initial values

    
end
