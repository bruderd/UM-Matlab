function s = state_encode(x, u, dxdt, params)
% converts states and inputs ICs into decision variable s
    x0 = params.x0;
    u0 = params.u0;
    N = params.N;
    m = params.m;
    n = params.n;

    s = zeros((N+1)*n + (N+1)*m + N*n,1);
    for i = 0 : N
        s(i*n + 1 : i*n + n) = x(:,i+1);
    end
        
    for j = 0 : N
        s(i*n+n+j*m + 1 : i*n+n+j*m + m) = u(:,j+1);
    end
    
    for k = 0 : N-1
        s(i*n+n+j*m+m + 1 : i*n+n+j*m+m + n) = dxdt(:,k+1);
    end
        
%    s = s';
end