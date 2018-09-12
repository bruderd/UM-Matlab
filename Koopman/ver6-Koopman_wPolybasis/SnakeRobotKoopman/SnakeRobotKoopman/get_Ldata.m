function Ldata = get_Ldata( U, params )
%get_Ldata: Calculate the approximate infinitesimal generator for the
%Koopman operator U
%   U = expm(L Ts). This function solves for L

% % eliminate complex eigenvalues of U so that it can be inverted
% [V, D, W] = eig(U);
% Dreal = D;
% for i = 1 : size(Dreal,2)
%    if Dreal(i,i) == 0
%        Dreal(i,i) = 1e-12;
%    end
% end
% logU = V * logm(Dreal) * W';
% 
% Ldata = (1/params.Ts) * logU;
% Ldata = real(Ldata);    % ignore complex components


Ldata = (1/params.Ts) * logm(U);

end