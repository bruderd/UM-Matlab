function Ldata = get_Ldata( U, params )
%get_Ldata: Calculate the approximate infinitesimal generator for the
%Koopman operator U
%   U = expm(L Ts). This function solves for L

Ldata = (1/params.Ts) * logm(U);

end