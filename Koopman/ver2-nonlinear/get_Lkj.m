function [ Lkj ] = get_Lkj( k, j, params )
%get_Lkj: Builds the operator L_j^k = p_k * \frac{\partial}{\partial x_j} 
%   See Eq. 13 in Mauroy and Gonclaves, 2016

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