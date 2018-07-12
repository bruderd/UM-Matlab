function [ Lkj ] = get_Lkj( k, j, params )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[N, psi, n] = deal(params.N, params.psi, params.n);

% create basis vectors for R^n
e = eye(n);


Lkj = zeros(N, N);
for i = 1 : N
    for l = 1 : N
        if psi(:,i) == psi(:,k) + psi(:,l) - e(:,j)
            Lkj(i,l) = psi(j,l);
        else
            Lkj(i,l) = 0;
        end
    end
end

end

