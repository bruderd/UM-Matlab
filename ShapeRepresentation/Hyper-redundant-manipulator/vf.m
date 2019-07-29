function Alphadot = vf( t , Alpha , u , params )
%vf: Explicit dynamics for robot
%   Alpha = [ alpha ; alphadot ];
%   Alphadot = [ alphadot ; alphaddot ];

alpha = Alpha( 1 : params.Nlinks );
alphadot = Alpha( params.Nlinks+1 : end );

Dq = get_massMatrix( alpha );
nonInert = get_nonInert( alpha , alphadot , u );

% solve for acceleration terms
alphaddot = - Dq \ nonInert;
% alphaddot = - pinv(Dq) * nonInert;

% add in the damping/input
% dampNinput = get_dampNinput( alpha , alphadot , u );
% alphaddot = alphaddot + dampNinput;

% define output
Alphadot = [ alphadot ; alphaddot ];

end

