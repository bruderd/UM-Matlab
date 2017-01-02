function [x_index, u_index] = index(k, s, params)

    N = params.N;
    m = params.m;
    n = params.n;
    
    x_index = s(k*n + 1 : k*n + n, 1);
    u_index = s((N+1)*n + k*m + 1 : (N+1)*n + k*m + m, 1);
    
end