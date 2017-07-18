function thetaSS = find_thetaSS(u, params)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here

[B, N, L, I, rho] = deal(params.B, params.N, params.L, params.I, params.rho);
[u1, u2] = deal(u(1), u(2));

% Without extraaneous solution
if u1 ~= u2
    thetaSS = ( 3*L*(u2+u1) - sqrt(3)*sqrt( B^2*(u2-u1)^2 + 12*L^2*u2*u1 ) ) / (3*rho*(u2-u1));
else
    thetaSS = 0;
end

end

