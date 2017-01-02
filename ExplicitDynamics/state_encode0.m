function s = state_encode0(params)
% converts states and inputs ICs into decision variable s
    x0 = params.x0;
    u0 = params.u0;
    N = params.N;
    m = params.m;
    n = params.n;

    s = [x0; rand(N*n,1); u0; rand(N*m,1)];  %intial decision variable s (includes x and u for all time steps). Initialized to all ones other than initial values
end
