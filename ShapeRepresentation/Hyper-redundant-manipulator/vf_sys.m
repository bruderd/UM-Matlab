function Alphadot = vf_sys( t , Alpha , u , sys )
%vf: Explicit dynamics for robot
%   Alpha = [ alpha ; alphadot ];
%   Alphadot = [ alphadot ; alphaddot ];

params = sys.params;

alpha = Alpha( 1 : params.Nlinks );
alphadot = Alpha( params.Nlinks+1 : end );

Dq = sys.get_massMatrix( alpha );
nonInert = sys.get_nonInert( alpha , alphadot , u );

% solve for acceleration terms
alphaddot = - Dq \ nonInert;
% alphaddot = - pinv(Dq) * nonInert;

% add in the damping/input
% dampNinput = get_dampNinput( alpha , alphadot , u );
% alphaddot = alphaddot + dampNinput;

% define output
Alphadot = [ alphadot ; alphaddot ];

end

