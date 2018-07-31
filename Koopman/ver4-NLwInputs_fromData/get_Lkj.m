function [ Lkj ] = get_Lkj( k, j, params )
%get_Lkj: Builds the operator L_j^k = p_k * \frac{\partial}{\partial x_j} 
%   See Eq. 13 in Mauroy and Gonclaves, 2016

[N, psi, naug] = deal(params.N, params.psi, params.naug);

% create basis vectors for R^n
e = eye(naug);


Lkj = zeros(N, params.N1);
for i = 1 : N
    for l = 1 : params.N1
        if psi(:,i) == psi(:,k) + psi(:,l) - e(:,j)
            Lkj(i,l) = psi(j,l);
        else
            Lkj(i,l) = 0;
        end
    end
end

end