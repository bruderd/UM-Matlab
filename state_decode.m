function [x, u] = state_decode(s, params)
% converts decision variable, s, to states and inputs
    N = params.N;
    m = params.m;
    n = params.n;
    
    x = zeros(N+1,n);
    u = zeros(N+1,m);
    for k = 1:(N + 1)
        x(k,:) = s( ( k - 1 )*n+1:k*n);
        u(k,:) = s(n*(N+1)+( k - 1 )*m+1:n*(N+1)+k*m);
    end
    
end