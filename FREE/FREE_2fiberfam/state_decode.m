function [x, u, dxdt] = state_decode(s, params)
% converts decision variable, s, to states and inputs
    N = params.N;
    m = params.m;
    n = params.n;
    
    x = zeros(N+1,n);
    u = zeros(N+1,m);
    dxdt = zeros(N,n);
    for k = 1:(N + 1)
        x(k,:) = s( (k-1)*n+1 : k*n);
        u(k,:) = s(n*(N+1)+(k-1)*m+1 : n*(N+1)+k*m);
    end
    
    for j = 1:N
        dxdt(j,:) = s((N+1)*n+(N+1)*m + (j-1)*n + 1 : (N+1)*n+(N+1)*m + (j-1)*n + n);
    end
    
    x = x';
    u = u';
    dxdt = dxdt';
end